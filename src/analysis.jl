module Analysis

# Orchestration module for post-pipeline analysis.
#
# Three analysis levels:
#   analyse_run    - per-run: pipeline summary CSV, taxa bar chart,
#                    filter composition, alpha diversity (PDFs)
#   analyse_group  - per-group: multi-run comparison charts, NMDS
#   analyse_study  - study-wide: cross-group NMDS, alpha comparison,
#                    PERMANOVA (when metadata.csv is present)
#
# analyse_run is called inside the @threads loop (no R calls).
# analyse_study is called once after the loop and dispatches to
# analyse_group internally (may call R for NMDS/PERMANOVA).
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export analyse_run, analyse_study, load_metadata

    using CSV, DataFrames, Dates, Logging, Statistics, YAML, RCall
    using ..PipelineTypes, ..PipelineLog, ..Config, ..DiversityMetrics, ..PipelinePlots

    ## Column identification
    """
        _sample_cols(df, db_meta) -> Vector{String}

    Return column names that hold per-sample ASV counts by excluding known
    taxonomy / metadata columns and any column ending with `_dada2` or
    `_boot`.
    """
    function _sample_cols(df::DataFrame, db_meta::DatabaseMeta)
        filter(names(df)) do col
            col ∉ db_meta.noncounts &&
            !endswith(col, "_dada2") &&
            !endswith(col, "_boot") &&
            !endswith(col, "_vsearch") &&
            !isempty(col)
        end
    end

    ##  CSV cache
    # Read a CSV once and reuse.  Stores (DataFrame, sample_cols).
    const _CSVCache = Dict{String, Tuple{DataFrame, Vector{String}}}

    function _cached_read(cache::_CSVCache, path::String, db_meta::DatabaseMeta)
        haskey(cache, path) && return cache[path]
        df    = CSV.read(path, DataFrame)
        scols = _sample_cols(df, db_meta)
        # Ensure sample columns are numeric (DADA2 tax_counts.csv may write them as strings).
        for col in scols
            if eltype(df[!, col]) <: AbstractString || eltype(df[!, col]) == Union{Missing, String}
                df[!, col] = [ismissing(v) ? missing : parse(Float64, v) for v in df[!, col]]
            elseif !(eltype(df[!, col]) <: Union{Missing, Number})
                df[!, col] = passmissing(x -> Float64(x)).(df[!, col])
            end
        end
        cache[path] = (df, scols)
        return df, scols
    end

    # Taxonomy methods: vsearch uses standard column names, dada2 uses _dada2 suffix.
    const _TAX_METHODS = ["vsearch", "dada2"]

    """
        _dada2_taxonomy_view(df, levels) -> DataFrame

    Return a copy where `_dada2` taxonomy columns are renamed to standard names
    (vsearch columns get `_vsearch` suffix). Allows plotting/reporting code to
    work unchanged on DADA2 taxonomy.
    """
    function _dada2_taxonomy_view(df::DataFrame, levels::Vector{String})
        vdf = copy(df)
        for rank in levels
            dada2_col = rank * "_dada2"
            dada2_col in names(vdf) || continue
            if rank in names(vdf)
                rename!(vdf, rank => rank * "_vsearch")
            end
            rename!(vdf, dada2_col => rank)
        end
        return vdf
    end

    # Check whether a DataFrame has any DADA2 taxonomy columns.
    _has_dada2(df::DataFrame, levels::Vector{String}) = any(c -> endswith(c, "_dada2") &&
        any(r -> startswith(c, r), levels), names(df))

    # For a given method, return the appropriate source keys from a MergedTables.
    # vsearch: keys without _dada2 suffix (existing behaviour)
    # dada2: for each filter stem, use "{stem}_dada2" key; for "merged" use "merged" (same CSV, both taxonomies)
    function _method_source_keys(merged::MergedTables, method::String)
        if method == "dada2"
            keys_out = String[]
            for k in sort(collect(keys(merged.tables)))
                k == "merged" && continue
                endswith(k, "_dada2") && push!(keys_out, k)
            end
            push!(keys_out, "merged")
            return keys_out
        else
            return sort([k for k in keys(merged.tables)
                         if k != "merged" && !endswith(k, "_dada2")]) |>
                   ks -> vcat(ks, ["merged"])
        end
    end

    # Get the CSV path and apply taxonomy view for a method.
    function _method_df(raw_df::DataFrame, method::String, levels::Vector{String})
        method == "dada2" ? _dada2_taxonomy_view(raw_df, levels) : raw_df
    end

    # Source dirname for method-qualified keys (strip _dada2 suffix for directory names).
    _method_source_dirname(key::String) = _source_dirname(
        endswith(key, "_dada2") ? key[1:end-6] : key)

    # Lowest assigned taxonomic rank for a row.
    function _lowest_rank_label(row, levels::Vector{String})
        label = "Unclassified"
        for rank in levels
            if hasproperty(row, Symbol(rank))
                val = row[Symbol(rank)]
                if !ismissing(val) && !isempty(strip(string(val)))
                    label = string(val)
                end
            end
        end
        return label
    end

    ## Stem-to-column matching
    # Match pipeline_stats sample stems to full FASTQ column names.
    # Returns a Dict mapping each stem to the matching column name
    # (or the stem itself if no match is found).
    function _stem_to_colname(stems::AbstractVector, col_names::Vector{String})
        mapping = Dict{String, String}()
        for stem in stems
            s = String(stem)
            matched = false
            for col in col_names
                if startswith(col, s * "_") || col == s
                    mapping[s] = col
                    matched = true
                    break
                end
            end
            matched || (mapping[s] = s)
        end
        return mapping
    end

    ## Source label helper
    _source_label(key::String) = key == "merged" ? "unfiltered" : key

    ## Taxa rank iteration
    # Derive taxa chart ranks from database levels and optional config override.
    function _taxa_ranks(levels::Vector{String}, analysis_cfg::Dict=Dict())
        configured = get(get(analysis_cfg, "taxa_bar", Dict()), "ranks", nothing)
        if !isnothing(configured) && configured isa Vector
            return Tuple{String, Union{String,Nothing}}[
                (lowercase(string(r)), string(r) == "asv" ? nothing : string(r)) for r in configured
            ]
        end
        n = length(levels)
        ranks = Tuple{String, Union{String,Nothing}}[("asv", nothing)]
        n >= 1 && push!(ranks, (lowercase(levels[end]),   levels[end]))
        n >= 2 && push!(ranks, (lowercase(levels[end-1]), levels[end-1]))
        return ranks
    end

    # Derive report ranks from database levels and optional config override.
    function _report_ranks(levels::Vector{String}, analysis_cfg::Dict=Dict())
        configured = get(get(analysis_cfg, "taxa_bar", Dict()), "report_ranks", nothing)
        if !isnothing(configured) && configured isa Vector
            return [(string(r), string(r)) for r in configured]
        end
        n = length(levels)
        ranks = Tuple{String,String}[]
        n >= 3 && push!(ranks, (levels[end-2], levels[end-2]))
        n >= 2 && push!(ranks, (levels[end-1], levels[end-1]))
        n >= 1 && push!(ranks, (levels[end],   levels[end]))
        return ranks
    end

    # Source key -> directory name for figures.
    _source_dirname(key::String) = key == "merged" ? "unfiltered" : key

    function _generate_taxa_charts(df::DataFrame, scols::Vector{String},
                                    figures_dir::String, top_n::Int,
                                    subtitle::Union{Nothing,String},
                                    source_key::String;
                                    ranks, rank_order::Vector{String})
        src_dir = _source_dirname(source_key)
        for (rankdir, rank) in ranks
            dir = joinpath(figures_dir, src_dir, rankdir)
            mkpath(dir)
            taxa_bar_chart(df, scols, joinpath(dir, "taxa_bar.pdf");
                           top_n, rank, relative=true, subtitle, rank_order)
            taxa_bar_chart(df, scols, joinpath(dir, "taxa_bar_absolute.pdf");
                           top_n, rank, relative=false, subtitle, rank_order)
        end
    end

    ## Top-taxa report section
    function _top_taxa_section(df::DataFrame, scols::Vector{String},
                                rank_name::String, rank_col::String; n::Int=20)
        sym = Symbol(rank_col)
        hasproperty(df, sym) || return ""

        # Label each row, sum counts across all samples.
        labels = String[]
        totals = Float64[]
        for row in eachrow(df)
            val = row[sym]
            label = (ismissing(val) || isempty(strip(string(val)))) ? "Unclassified" : string(val)
            push!(labels, label)
            push!(totals, sum(col -> begin
                v = row[Symbol(col)]
                ismissing(v) ? 0.0 : Float64(v)
            end, scols))
        end

        # Aggregate by label.
        agg = Dict{String, Float64}()
        for (l, t) in zip(labels, totals)
            agg[l] = get(agg, l, 0.0) + t
        end

        sorted = sort(collect(agg); by=last, rev=true)
        grand_total = sum(last, sorted; init=0.0)

        buf = IOBuffer()
        print(buf, rpad("Rank", 4), rpad(rank_name, 30), rpad("Reads", 14), "Percent\n")
        for (i, (label, count)) in enumerate(sorted)
            i > n && break
            pct = grand_total > 0 ? round(100.0 * count / grand_total; digits=1) : 0.0
            print(buf, rpad(string(i), 4), rpad(label, 30),
                  rpad(string(Int(count)), 14), "$(pct)%\n")
        end
        return String(take!(buf))
    end

    ## Report generator (dual: filtered + merged)
    function _generate_report(df::DataFrame, scols::Vector{String},
                               stats_df, merged_df::DataFrame,
                               src_key::String, report_path::String,
                               run_name::String;
                               stats_key::String=src_key,
                               report_ranks::Vector{Tuple{String,String}}=[("Family","Family"),("Genus","Genus"),("Species","Species")])
        src_label = _source_label(src_key)
        report_sections = Pair{String, String}[]

        # Pipeline statistics table.
        if !isnothing(stats_df)
            buf = IOBuffer()
            _print_stats_table(buf, stats_df)
            if "input" in names(stats_df)
                total_input = sum(stats_df.input)
                # Use the reads column matching this report's stats key.
                target_col = "reads_" * stats_key
                if target_col in names(stats_df)
                    total_final = sum(stats_df[!, target_col])
                    pct = total_input > 0 ? round(100.0 * total_final / total_input; digits=1) : 0.0
                    println(buf, "\nOverall retention (input -> $(src_label)): $(pct)%")
                end
            end
            push!(report_sections, "Pipeline Statistics" => String(take!(buf)))
        end

        # ASV summary.
        merged_asvs = nrow(merged_df)
        filter_asvs = nrow(df)
        buf = IOBuffer()
        println(buf, "Total ASVs (merged):    $merged_asvs")
        if src_key != "merged"
            println(buf, "ASVs after filter:      $filter_asvs")
        end
        push!(report_sections, "ASV Summary" => String(take!(buf)))

        # Alpha diversity table.
        alpha = _compute_alpha(df, scols)
        buf = IOBuffer()
        _print_alpha_table(buf, alpha)
        push!(report_sections, "Alpha Diversity" => String(take!(buf)))

        # Top-20 taxa tables.
        for (rank_name, rank_col) in report_ranks
            section = _top_taxa_section(df, scols, rank_name, rank_col; n=20)
            !isempty(section) && push!(report_sections, "Top 20 $(rank_name) ($(src_label))" => section)
        end

        _write_report(report_path,
                      "Analysis Report: $run_name\n  Source: $src_label",
                      report_sections)
    end

    ## Text report writer
    function _write_report(path::String, title::String,
                           sections::Vector{Pair{String, String}})
        open(path, "w") do io
            sep = "=" ^ 64
            println(io, sep)
            println(io, "  ", title)
            println(io, "  Generated: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
            println(io, sep)
            for (heading, body) in sections
                println(io)
                println(io, "--- ", heading, " ---")
                println(io, body)
            end
        end
        @info "Written: $path"
    end

    ## Source table selection
    # Pick the primary analysis source from a MergedTables: the first
    # filter CSV (alphabetically) if any filters exist, else "merged".
    function _source_key(merged::MergedTables)
        "merged"
    end

    ## Priority-based filter composition
    # Assign each ASV to the first matching filter (by pipeline.yml order).
    # Returns (filter_totals, sample_cols) or (nothing, nothing).
    function _priority_filter_composition(merged::MergedTables, cache::_CSVCache,
                                           db_meta::DatabaseMeta)
        merged_csv = merged.tables["merged"]
        isfile(merged_csv) || return nothing, nothing
        merged_df, merged_scols = _cached_read(cache, merged_csv, db_meta)
        isempty(merged_scols) && return nothing, nothing

        # Use filter_order for priority; fall back to sorted keys.
        filter_keys = !isempty(merged.filter_order) ? merged.filter_order :
                      sort([k for k in keys(merged.tables) if k != "merged"])
        isempty(filter_keys) && return nothing, nothing

        # Per-sample totals from merged (unfiltered).
        merged_totals = Float64[
            sum(v -> ismissing(v) ? 0.0 : Float64(v), merged_df[!, s])
            for s in merged_scols
        ]

        # Load filter DataFrames.
        filter_dfs = Dict{String, DataFrame}()
        for fk in filter_keys
            haskey(merged.tables, fk) || continue
            fpath = merged.tables[fk]
            (isfile(fpath) && filesize(fpath) > 0) || continue
            filter_dfs[fk] = first(_cached_read(cache, fpath, db_meta))
        end

        # Track which SeqNames have been claimed.
        claimed = Set{String}()
        filter_totals = Dict{String, Vector{Float64}}()

        for fk in filter_keys
            haskey(filter_dfs, fk) || continue
            fdf = filter_dfs[fk]
            totals = zeros(Float64, length(merged_scols))
            for row in eachrow(fdf)
                seq = string(row.SeqName)
                seq in claimed && continue
                push!(claimed, seq)
                for (i, s) in enumerate(merged_scols)
                    v = hasproperty(row, Symbol(s)) ? row[Symbol(s)] : missing
                    totals[i] += ismissing(v) ? 0.0 : Float64(v)
                end
            end
            filter_totals[fk] = totals
        end

        # Unclassified = merged totals minus all attributed.
        attributed = zeros(Float64, length(merged_scols))
        for vals in values(filter_totals)
            attributed .+= vals
        end
        unclassified = max.(merged_totals .- attributed, 0.0)
        if any(>(0), unclassified)
            filter_totals["Unclassified"] = unclassified
        end

        return filter_totals, merged_scols
    end

    ## Count-matrix builders
    # Extract a samples x features count matrix from a single DataFrame.
    # Rows = samples (one per element of `scols`), columns = ASV rows.
    function _counts_matrix(df::DataFrame, scols::Vector{String})
        n_samples  = length(scols)
        n_features = nrow(df)
        mat = zeros(Float64, n_samples, n_features)
        for (j, row) in enumerate(eachrow(df))
            for (i, col) in enumerate(scols)
                v = row[Symbol(col)]
                mat[i, j] = ismissing(v) ? 0.0 : Float64(v)
            end
        end
        return mat
    end

    # Build a taxonomy-aggregated count matrix across multiple runs.
    #
    # ASV identifiers (seq1, seq2, ...) are local to each run and cannot
    # be compared directly.  This function aggregates counts to the
    # lowest assigned taxonomic rank, producing a shared feature space
    # suitable for between-run Bray-Curtis / NMDS.
    #
    # Returns (matrix, all_sample_names, taxon_labels).
    function _build_combined_counts(
        dfs::Vector{DataFrame},
        scols_per_df::Vector{Vector{String}},
        levels::Vector{String},
    )
        # taxon -> sample -> accumulated count
        taxa_counts = Dict{String, Dict{String, Float64}}()
        all_samples = String[]

        for (df, scols) in zip(dfs, scols_per_df)
            append!(all_samples, scols)
            for row in eachrow(df)
                label = _lowest_rank_label(row, levels)
                td = get!(taxa_counts, label, Dict{String, Float64}())
                for col in scols
                    v = row[Symbol(col)]
                    td[col] = get(td, col, 0.0) + (ismissing(v) ? 0.0 : Float64(v isa AbstractString ? parse(Float64, v) : v))
                end
            end
        end

        taxa_labels = sort(collect(keys(taxa_counts)))
        n_samples   = length(all_samples)
        n_taxa      = length(taxa_labels)
        mat = zeros(Float64, n_samples, n_taxa)

        for (j, taxon) in enumerate(taxa_labels)
            td = taxa_counts[taxon]
            for (i, sample) in enumerate(all_samples)
                mat[i, j] = get(td, sample, 0.0)
            end
        end

        return mat, all_samples, taxa_labels
    end

    ## Alpha diversity helper
    # Compute per-sample alpha diversity from a merged/filtered CSV.
    function _compute_alpha(df::DataFrame, scols::Vector{String})
        out = DataFrame(sample=String[], richness=Int[],
                        shannon=Float64[], simpson=Float64[])
        for col in scols
            counts = [ismissing(v) ? 0 : Int(v isa AbstractString ? parse(Int, v) : round(Int, Float64(v))) for v in df[!, col]]
            push!(out, (col, richness(counts), shannon(counts), simpson(counts)))
        end
        return out
    end

    # Metadata
    """
        load_metadata(start_dir, study_dir) -> Union{DataFrame, Nothing}

    Walk from `start_dir` upward to `study_dir` (inclusive), returning the
    first `metadata.csv` found as a DataFrame.  Returns `nothing` when no
    metadata file exists at any level.
    """
    function load_metadata(start_dir::String, study_dir::String)
        dir  = abspath(start_dir)
        stop = abspath(study_dir)
        while true
            csv = joinpath(dir, "metadata.csv")
            isfile(csv) && return CSV.read(csv, DataFrame)
            dir == stop && break
            parent = dirname(dir)
            parent == dir && break          # filesystem root
            dir = parent
        end
        return nothing
    end

    # R: NMDS + PERMANOVA
    # NMDS via vegan::metaMDS.
    # `mat` is samples x features (community matrix).
    # Returns (coords::Matrix{Float64}[nx2], stress::Float64).
    # On failure returns a NaN-filled matrix and NaN stress.
    function _run_nmds(mat::Matrix{Float64}, r_lock::ReentrantLock)
        lock(r_lock) do
            @rput mat
            R"""
            suppressPackageStartupMessages(library(vegan))
            set.seed(42)
            nmds_res <- tryCatch(
                metaMDS(mat, distance = "bray", k = 2, trymax = 200,
                        autotransform = FALSE, trace = 0),
                error = function(e) NULL
            )
            if (!is.null(nmds_res)) {
                nmds_coords <- nmds_res$points
                nmds_stress <- nmds_res$stress
            } else {
                nmds_coords <- matrix(NA_real_, nrow = nrow(mat), ncol = 2)
                nmds_stress <- NA_real_
            }
            """
            coords = rcopy(R"nmds_coords")::Matrix{Float64}
            stress = rcopy(R"nmds_stress")::Float64
            return coords, stress
        end
    end

    # PERMANOVA via vegan::adonis2.
    # Returns the captured text output, or `nothing` on failure.
    function _run_permanova(mat::Matrix{Float64}, metadata::DataFrame,
                            r_lock::ReentrantLock)
        covariates = [c for c in names(metadata) if lowercase(c) != "sample"]
        isempty(covariates) && return nothing
        formula_rhs = join(covariates, " + ")

        lock(r_lock) do
            meta_r = copy(metadata)
            @rput mat meta_r formula_rhs
            R"""
            suppressPackageStartupMessages(library(vegan))
            set.seed(42)
            dist_mat <- vegdist(mat, method = "bray")
            form     <- as.formula(paste("dist_mat ~", formula_rhs))
            perm_res <- tryCatch(
                adonis2(form, data = meta_r, permutations = 999),
                error = function(e) NULL
            )
            if (!is.null(perm_res)) {
                perm_text <- paste(capture.output(print(perm_res)), collapse = "\n")
            } else {
                perm_text <- NA_character_
            }
            """
            txt = rcopy(R"perm_text")
            return (ismissing(txt) || txt == "NA") ? nothing : txt
        end
    end

    ## Pipeline summary CSV
    function _pipeline_summary(project::ProjectCtx, merged::MergedTables,
                               db_meta::DatabaseMeta)
        analysis_dir = joinpath(project.dir, "analysis")
        summary_path = joinpath(analysis_dir, "pipeline_summary.csv")

        # Read DADA2's pipeline_stats.csv.
        # R's write.csv writes row names as the first (unnamed) column.
        stats_csv = joinpath(project.dir, "dada2", "Tables", "pipeline_stats.csv")
        if !isfile(stats_csv)
            @warn "pipeline_stats.csv not found at $stats_csv - skipping pipeline_summary"
            return nothing
        end

        stats = CSV.read(stats_csv, DataFrame)
        first_col = names(stats)[1]
        if first_col != "sample"
            rename!(stats, first_col => "sample")
        end

        # Add per-sample read totals from each merged table.
        # Stems in pipeline_stats (e.g. "JIN-Nu-ves55") must be matched to
        # full FASTQ column names (e.g. "JIN-Nu-ves55_R1_filt.fastq.gz").
        for key in sort(collect(keys(merged.tables)))
            csv_path = merged.tables[key]
            isfile(csv_path) || continue
            df    = CSV.read(csv_path, DataFrame)
            scols = _sample_cols(df, db_meta)

            totals = Dict{String, Int}()
            for col in scols
                totals[col] = sum(v -> ismissing(v) ? 0 : Int(v), df[!, col]; init=0)
            end

            # Match stems to full column names.
            stem_map = _stem_to_colname(stats.sample, scols)

            # Prefix with "reads_" to avoid colliding with DADA2's own columns.
            col_name = "reads_" * key
            stats[!, col_name] = [get(totals, get(stem_map, String(s), ""), 0)
                                  for s in stats.sample]
        end

        mkpath(analysis_dir)
        CSV.write(summary_path, stats)
        @info "Written: $summary_path"
        log_written(project, summary_path)
        return stats
    end

    ## Level 1 - Per-run analysis
    """
        analyse_run(project, merged, asvs, db_meta; plot_lock=nothing)

    Produce per-run analysis outputs under `{project.dir}/analysis/`.
    When `plot_lock` is provided, CairoMakie calls are serialized behind
    it so multiple runs can prepare data in parallel.
    """
    function analyse_run(project::ProjectCtx, merged::MergedTables,
                         asvs::ASVResult, db_meta::DatabaseMeta;
                         plot_lock::Union{Nothing,ReentrantLock}=nothing)
        analysis_dir = joinpath(project.dir, "analysis")
        figures_dir  = joinpath(analysis_dir, "Figures")
        merged_csv   = merged.tables["merged"]

        filter_keys  = sort([k for k in keys(merged.tables) if k != "merged"])
        has_filters  = !isempty(filter_keys)
        all_source_keys = vcat(filter_keys, ["merged"])

        # Output paths.
        summary_path  = joinpath(analysis_dir, "pipeline_summary.csv")
        stages_pdf    = joinpath(figures_dir,  "pipeline_stages.pdf")
        alpha_pdf     = joinpath(figures_dir,  "alpha_diversity.pdf")
        filter_pdf    = joinpath(figures_dir,  "filter_composition.pdf")

        # Read analysis config from the cascade (needed for skip guard and charts).
        config_path  = write_run_config(project)
        analysis_cfg = get(YAML.load_file(config_path), "analysis", Dict())
        taxa_cfg     = get(analysis_cfg, "taxa_bar", Dict())
        top_n        = get(taxa_cfg, "top_n", 15)
        taxa_ranks   = _taxa_ranks(db_meta.levels, analysis_cfg)
        report_ranks = _report_ranks(db_meta.levels, analysis_cfg)

        required_outputs = String[summary_path, stages_pdf, alpha_pdf]
        has_filters && push!(required_outputs, filter_pdf)
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(merged, method)
            for src in msrc_keys
                sd = _method_source_dirname(src)
                for (rankdir, _) in taxa_ranks
                    push!(required_outputs, joinpath(figures_dir, method, sd, rankdir, "taxa_bar.pdf"))
                    push!(required_outputs, joinpath(figures_dir, method, sd, rankdir, "taxa_bar_absolute.pdf"))
                end
            end
            for src in msrc_keys
                push!(required_outputs, joinpath(analysis_dir, "analysis_report_$(method)_$(_method_source_dirname(src)).txt"))
            end
        end

        # skip guard
        merged_mtime = isfile(merged_csv) ? mtime(merged_csv) : time()
        if all(isfile, required_outputs) &&
           all(f -> mtime(f) > merged_mtime, required_outputs)
            @info "Skipping analyse_run: outputs up to date in $analysis_dir"
            return
        end

        mkpath(figures_dir)

        # Helper: serialize CairoMakie calls if plot_lock is provided.
        _plot(f) = isnothing(plot_lock) ? f() : lock(f, plot_lock)

        # CSV cache
        cache = _CSVCache()

        src_key   = _source_key(merged)
        src_label = _source_label(src_key)
        subtitle  = "Source: $src_label"

        # pipeline summary
        stats_df = _pipeline_summary(project, merged, db_meta)

        # pipeline stages plot
        if !isnothing(stats_df)
            _plot() do
                pipeline_stats_plot(stats_df, stages_pdf; subtitle)
            end
            @info "Written: $stages_pdf"
            log_written(project, stages_pdf)
        end

        # taxa bar charts for all sources at multiple ranks (dual method)
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(merged, method)
            for src in msrc_keys
                # For dada2 method with "merged" key, use merged.csv;
                # for _dada2 suffixed keys, use the corresponding CSV.
                csv_key = src
                src_csv_path = get(merged.tables, csv_key, "")
                if isempty(src_csv_path) || !isfile(src_csv_path) || filesize(src_csv_path) == 0
                    # For "merged" key under dada2 method, fall back to merged.csv
                    src == "merged" || continue
                    src_csv_path = merged.tables["merged"]
                    (!isfile(src_csv_path) || filesize(src_csv_path) == 0) && continue
                end
                raw_df, src_scols = _cached_read(cache, src_csv_path, db_meta)
                isempty(src_scols) && continue
                # Skip dada2 method if no dada2 taxonomy columns
                method == "dada2" && !_has_dada2(raw_df, db_meta.levels) && continue
                view_df = _method_df(raw_df, method, db_meta.levels)
                sd = _method_source_dirname(src)
                src_sub = "Source: $sd ($method)"
                method_fig_dir = joinpath(figures_dir, method)
                _plot() do
                    _generate_taxa_charts(view_df, src_scols, method_fig_dir,
                                           top_n, src_sub, src == "merged" ? "merged" :
                                           (endswith(src, "_dada2") ? src[1:end-6] : src);
                                           ranks=taxa_ranks, rank_order=db_meta.levels)
                end
            end
        end

        # filter composition (priority-based)
        if has_filters
            filter_totals, comp_scols = _priority_filter_composition(merged, cache, db_meta)
            if !isnothing(filter_totals) && !isempty(filter_totals)
                _plot() do
                    filter_composition_plot(filter_totals, comp_scols, filter_pdf;
                                            subtitle, colour_overrides=merged.filter_colours)
                end
                @info "Written: $filter_pdf"
                log_written(project, filter_pdf)
            end
        end

        # alpha diversity (from primary source)
        primary_csv = merged.tables[src_key]
        if isfile(primary_csv) && filesize(primary_csv) > 0
            prim_df, prim_scols = _cached_read(cache, primary_csv, db_meta)
            if !isempty(prim_scols)
                alpha = _compute_alpha(prim_df, prim_scols)
                _plot() do
                    alpha_diversity_plot(alpha, alpha_pdf; subtitle)
                end
                @info "Written: $alpha_pdf"
                log_written(project, alpha_pdf)
            end
        end

        # analysis reports (dual method, one per source per method)
        run_name = basename(project.dir)
        merged_df_full, _ = _cached_read(cache, merged_csv, db_meta)
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(merged, method)
            for src in msrc_keys
                csv_key = src
                src_csv_path = get(merged.tables, csv_key, "")
                if isempty(src_csv_path) || !isfile(src_csv_path) || filesize(src_csv_path) == 0
                    src == "merged" || continue
                    src_csv_path = merged.tables["merged"]
                    (!isfile(src_csv_path) || filesize(src_csv_path) == 0) && continue
                end
                raw_df, src_scols = _cached_read(cache, src_csv_path, db_meta)
                isempty(src_scols) && continue
                method == "dada2" && !_has_dada2(raw_df, db_meta.levels) && continue
                view_df = _method_df(raw_df, method, db_meta.levels)
                sd = _method_source_dirname(src)
                report_path = joinpath(analysis_dir, "analysis_report_$(method)_$(sd).txt")
                clean_src = src == "merged" ? "merged" : (endswith(src, "_dada2") ? src[1:end-6] : src)
                _generate_report(view_df, src_scols, stats_df, merged_df_full,
                                 clean_src, report_path, run_name;
                                 stats_key=src, report_ranks=report_ranks)
                log_written(project, report_path)
            end
        end
    end

    # Format a pipeline stats DataFrame as a text table.
    function _print_stats_table(io::IO, stats_df::DataFrame)
        cols = names(stats_df)
        # Header.
        print(io, rpad("Sample", 20))
        for col in cols
            col == "sample" && continue
            print(io, rpad(col, 16))
        end
        println(io)
        # Rows.
        for row in eachrow(stats_df)
            print(io, rpad(PipelinePlots._display_name(String(row.sample)), 20))
            for col in cols
                col == "sample" && continue
                v = row[Symbol(col)]
                print(io, rpad(ismissing(v) ? "0" : string(Int(v)), 16))
            end
            println(io)
        end
    end

    # Format an alpha diversity DataFrame as a text table.
    function _print_alpha_table(io::IO, alpha_df::DataFrame)
        print(io, rpad("Sample", 20))
        println(io, rpad("Richness", 12), rpad("Shannon", 12), "Simpson")
        for row in eachrow(alpha_df)
            print(io, rpad(PipelinePlots._display_name(String(row.sample)), 20))
            print(io, rpad(string(row.richness), 12))
            print(io, rpad(string(round(row.shannon; digits=3)), 12))
            println(io, round(row.simpson; digits=3))
        end
    end

    ## Level 2 - Per-group analysis
    function _analyse_group(group_dir::String,
                            members::Vector{Tuple{ProjectCtx, MergedTables}},
                            r_lock::ReentrantLock,
                            db_meta::DatabaseMeta;
                            plot_lock::Union{Nothing,ReentrantLock}=nothing)
        analysis_dir = joinpath(group_dir, "analysis")
        figures_dir  = joinpath(analysis_dir, "Figures")

        alpha_pdf     = joinpath(figures_dir, "alpha_comparison.pdf")
        nmds_pdf      = joinpath(figures_dir, "nmds.pdf")
        filter_pdf    = joinpath(figures_dir, "filter_composition.pdf")

        # Helper: serialize CairoMakie calls if plot_lock is provided.
        _plot(f) = isnothing(plot_lock) ? f() : lock(f, plot_lock)

        # Determine source keys.
        src_key   = _source_key(members[1][2])
        src_label = _source_label(src_key)
        subtitle  = "Source: $src_label"

        # Derive taxa/report ranks from db_meta.
        analysis_cfg = Dict()
        try
            config_path = write_run_config(members[1][1])
            analysis_cfg = get(YAML.load_file(config_path), "analysis", Dict())
        catch; end
        taxa_ranks   = _taxa_ranks(db_meta.levels, analysis_cfg)
        report_ranks = _report_ranks(db_meta.levels, analysis_cfg)

        filter_keys = sort([k for k in keys(members[1][2].tables) if k != "merged"])
        has_any_filters = any(pair -> length(pair[2].tables) > 1, members)
        all_source_keys = vcat(filter_keys, ["merged"])

        # skip guard
        newest_merged = maximum(
            mtime(m.tables["merged"])
            for (_, m) in members if isfile(m.tables["merged"]);
            init=0.0
        )
        required_outputs = String[alpha_pdf, nmds_pdf]
        has_any_filters && push!(required_outputs, filter_pdf)
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(members[1][2], method)
            for src in msrc_keys
                sd = _method_source_dirname(src)
                for (rankdir, _) in taxa_ranks
                    push!(required_outputs, joinpath(figures_dir, method, sd, rankdir, "taxa_bar.pdf"))
                    push!(required_outputs, joinpath(figures_dir, method, sd, rankdir, "taxa_bar_absolute.pdf"))
                    push!(required_outputs, joinpath(figures_dir, method, sd, rankdir, "group_comparison.pdf"))
                    push!(required_outputs, joinpath(figures_dir, method, sd, rankdir, "group_comparison_absolute.pdf"))
                end
            end
            for src in msrc_keys
                push!(required_outputs, joinpath(analysis_dir, "analysis_report_$(method)_$(_method_source_dirname(src)).txt"))
            end
        end

        if all(isfile, required_outputs) &&
           all(f -> mtime(f) > newest_merged, required_outputs)
            @info "Skipping analyse_group: outputs up to date in $figures_dir"
            return
        end

        reset_log(group_dir)
        mkpath(figures_dir)

        # CSV cache
        cache = _CSVCache()

        # helper: collect per-run data for a given source key
        function _collect_source_data(source::String)
            dfs    = DataFrame[]
            scols  = Vector{String}[]
            names_ = String[]
            labels = String[]
            all_s  = String[]
            for (proj, merged) in members
                haskey(merged.tables, source) || continue
                csv = merged.tables[source]
                isfile(csv) || continue
                df, sc = _cached_read(cache, csv, db_meta)
                isempty(sc) && continue
                rn = basename(proj.dir)
                push!(dfs, df); push!(scols, sc); push!(names_, rn)
                append!(all_s, sc); append!(labels, fill(rn, length(sc)))
            end
            return dfs, scols, names_, labels, all_s
        end

        # Primary source data (for NMDS, alpha, filter composition).
        run_dfs, run_scols, run_names, run_labels, all_scols =
            _collect_source_data(src_key)
        isempty(run_dfs) && return

        # taxa bar & group comparison for ALL sources (dual method)
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(members[1][2], method)
            for src in msrc_keys
                s_dfs, s_scols, s_names, _, s_all = _collect_source_data(
                    endswith(src, "_dada2") ? src : src)
                # For "merged" under dada2, use the same "merged" key
                if isempty(s_dfs) && method == "dada2" && src == "merged"
                    s_dfs, s_scols, s_names, _, s_all = _collect_source_data("merged")
                end
                isempty(s_dfs) && continue
                # Skip dada2 if no dada2 columns
                method == "dada2" && !_has_dada2(s_dfs[1], db_meta.levels) && continue
                view_dfs = [_method_df(df, method, db_meta.levels) for df in s_dfs]
                combined = reduce((a, b) -> vcat(a, b; cols=:union), view_dfs)
                sd = _method_source_dirname(src)
                src_sub = "Source: $sd ($method)"
                method_fig_dir = joinpath(figures_dir, method)
                # Use cleaned source key for directory naming
                clean_src = src == "merged" ? "merged" :
                            (endswith(src, "_dada2") ? src[1:end-6] : src)
                _plot() do
                    _generate_taxa_charts(combined, s_all, method_fig_dir, 15, src_sub, clean_src;
                                           ranks=taxa_ranks, rank_order=db_meta.levels)
                    for (rankdir, rank) in taxa_ranks
                        dir = joinpath(method_fig_dir, _source_dirname(clean_src), rankdir)
                        mkpath(dir)
                        group_comparison_chart(view_dfs, s_scols, s_names,
                                               joinpath(dir, "group_comparison.pdf");
                                               top_n=15, rank, relative=true,
                                               subtitle=src_sub,
                                               rank_order=db_meta.levels)
                        group_comparison_chart(view_dfs, s_scols, s_names,
                                               joinpath(dir, "group_comparison_absolute.pdf");
                                               top_n=15, rank, relative=false,
                                               subtitle=src_sub,
                                               rank_order=db_meta.levels)
                    end
                end
            end
        end

        # group alpha diversity (boxplot by run)
        all_alpha    = DataFrame[]
        alpha_labels = String[]
        for (df, scols, rname) in zip(run_dfs, run_scols, run_names)
            alpha = _compute_alpha(df, scols)
            push!(all_alpha, alpha)
            push!(alpha_labels, rname)
        end
        if !isempty(all_alpha)
            _plot() do
                alpha_boxplot(all_alpha, alpha_labels, alpha_pdf; subtitle)
            end
            @info "Written: $alpha_pdf"
            log_written(group_dir, alpha_pdf)
        end

        # group NMDS
        mat, _, _ = _build_combined_counts(run_dfs, run_scols, db_meta.levels)
        nmds_stress = NaN
        if size(mat, 1) >= 3
            coords, stress = _run_nmds(mat, r_lock)
            if !any(isnan, coords)
                nmds_cfg = Dict()
                try
                    config_path = write_run_config(members[1][1])
                    nmds_cfg = get(get(YAML.load_file(config_path), "analysis", Dict()),
                                   "nmds", Dict())
                catch; end
                max_stress = get(nmds_cfg, "max_stress", 0.2)
                if !isnan(stress) && stress > max_stress
                    @warn "NMDS stress $(round(stress; digits=3)) exceeds threshold $max_stress for group $(basename(group_dir))"
                    pipeline_log(group_dir, "WARN: NMDS stress $(round(stress; digits=3)) exceeds threshold $max_stress")
                end
                nmds_stress = stress
                _plot() do
                    nmds_plot(coords, all_scols, nmds_pdf;
                              colour_by=run_labels, stress=stress, subtitle)
                end
                @info "Written: $nmds_pdf"
                log_written(group_dir, nmds_pdf)
            else
                @warn "NMDS failed for group $(basename(group_dir))"
                pipeline_log(group_dir, "WARN: NMDS failed for group $(basename(group_dir))")
            end
        else
            @warn "Too few samples ($(size(mat, 1))) for NMDS in group $(basename(group_dir))"
        end

        # group filter composition
        if has_any_filters
            group_filter_totals  = Dict{String, Vector{Float64}}()
            group_sample_names   = String[]
            group_filter_colours = Dict{String,String}()
            for (proj, merged) in members
                ft, scols_m = _priority_filter_composition(merged, cache, db_meta)
                isnothing(ft) && continue
                append!(group_sample_names, scols_m)
                merge!(group_filter_colours, merged.filter_colours)
                for (fk, vals) in ft
                    existing = get(group_filter_totals, fk, Float64[])
                    group_filter_totals[fk] = vcat(existing, vals)
                end
            end

            if haskey(group_filter_totals, "Unclassified") &&
               !any(>(0), group_filter_totals["Unclassified"])
                delete!(group_filter_totals, "Unclassified")
            end

            if !isempty(group_filter_totals) && !isempty(group_sample_names)
                _plot() do
                    filter_composition_plot(group_filter_totals, group_sample_names, filter_pdf;
                                            subtitle, colour_overrides=group_filter_colours)
                end
                @info "Written: $filter_pdf"
                log_written(group_dir, filter_pdf)
            end
        end

        # analysis reports (dual method, per source, each with top-20 tables)
        # Shared sections: Per-Run Summary + NMDS.
        shared_sections = Pair{String, String}[]

        buf = IOBuffer()
        print(buf, rpad("Run", 20))
        println(buf, rpad("Samples", 10), rpad("Mean richness", 16),
                rpad("Mean Shannon", 16), "Mean Simpson")
        for (adf, rname) in zip(all_alpha, alpha_labels)
            n = nrow(adf)
            mr = round(mean(adf.richness); digits=0)
            ms = round(mean(adf.shannon); digits=3)
            mp = round(mean(adf.simpson); digits=3)
            if n > 1
                sr = round(std(adf.richness); digits=0)
                ss = round(std(adf.shannon); digits=2)
                sp = round(std(adf.simpson); digits=2)
                print(buf, rpad(rname, 20))
                println(buf, rpad(string(n), 10),
                        rpad("$(Int(mr)) +/- $(Int(sr))", 16),
                        rpad("$ms +/- $ss", 16),
                        "$mp +/- $sp")
            else
                print(buf, rpad(rname, 20))
                println(buf, rpad(string(n), 10),
                        rpad(string(Int(mr)), 16),
                        rpad(string(ms), 16),
                        string(mp))
            end
        end
        push!(shared_sections, "Per-Run Summary" => String(take!(buf)))

        if !isnan(nmds_stress)
            quality = nmds_stress < 0.2 ? "Good: below 0.2 threshold" : "Poor: above 0.2 threshold"
            push!(shared_sections, "NMDS" =>
                  "Stress: $(round(nmds_stress; digits=3)) ($quality)")
        end

        # Write one report per source per method.
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(members[1][2], method)
            for src in msrc_keys
                s_dfs, s_scols, _, _, s_all = _collect_source_data(
                    endswith(src, "_dada2") ? src : src)
                if isempty(s_dfs) && method == "dada2" && src == "merged"
                    s_dfs, s_scols, _, _, s_all = _collect_source_data("merged")
                end
                isempty(s_dfs) && continue
                method == "dada2" && !_has_dada2(s_dfs[1], db_meta.levels) && continue
                view_dfs = [_method_df(df, method, db_meta.levels) for df in s_dfs]
                combined = reduce((a, b) -> vcat(a, b; cols=:union), view_dfs)
                sd = _method_source_dirname(src)
                sl = sd == "unfiltered" ? "unfiltered" : sd

                sections = copy(shared_sections)
                for (rank_name, rank_col) in report_ranks
                    section = _top_taxa_section(combined, s_all, rank_name, rank_col; n=20)
                    !isempty(section) && push!(sections, "Top 20 $(rank_name) ($(sl))" => section)
                end

                rpath = joinpath(analysis_dir, "analysis_report_$(method)_$(sd).txt")
                _write_report(rpath,
                              "Comparison Report\n  Runs: $(join(run_names, ", "))\n  Source: $sl ($method)",
                              sections)
                log_written(group_dir, rpath)
            end
        end
    end

    ## Level 3 - Study-level analysis
    function _analyse_study_level(study_dir::String,
                                  valid::Vector{Tuple{ProjectCtx, MergedTables, DatabaseMeta}},
                                  r_lock::ReentrantLock,
                                  default_db_meta::DatabaseMeta;
                                  plot_lock::Union{Nothing,ReentrantLock}=nothing)
        analysis_dir = joinpath(study_dir, "analysis")
        figures_dir  = joinpath(analysis_dir, "Figures")
        nmds_pdf     = joinpath(figures_dir, "nmds.pdf")
        alpha_pdf    = joinpath(figures_dir, "alpha_comparison.pdf")
        perm_txt     = joinpath(analysis_dir, "permanova.txt")

        src_key   = _source_key(valid[1][2])
        filter_keys = sort([k for k in keys(valid[1][2].tables) if k != "merged"])

        # Derive taxa/report ranks from the first project's db_meta.
        analysis_cfg = Dict()
        try
            config_path = write_run_config(valid[1][1])
            analysis_cfg = get(YAML.load_file(config_path), "analysis", Dict())
        catch; end
        taxa_ranks   = _taxa_ranks(default_db_meta.levels, analysis_cfg)
        report_ranks = _report_ranks(default_db_meta.levels, analysis_cfg)

        # skip guard
        newest_merged = maximum(
            mtime(m.tables["merged"])
            for (_, m) in valid if isfile(m.tables["merged"]);
            init=0.0
        )
        required_outputs = [nmds_pdf, alpha_pdf]
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(valid[1][2], method)
            for src in msrc_keys
                push!(required_outputs, joinpath(analysis_dir, "analysis_report_$(method)_$(_method_source_dirname(src)).txt"))
            end
        end
        if all(isfile, required_outputs) &&
           all(f -> mtime(f) > newest_merged, required_outputs)
            @info "Skipping study-level analysis: outputs up to date"
            return
        end

        reset_log(study_dir)
        mkpath(figures_dir)

        # Helper: serialize CairoMakie calls if plot_lock is provided.
        _plot(f) = isnothing(plot_lock) ? f() : lock(f, plot_lock)

        cache = _CSVCache()

        src_label = _source_label(src_key)
        subtitle  = "Source: $src_label"

        # group projects by parent dir
        groups = Dict{String, Vector{Tuple{ProjectCtx, MergedTables, DatabaseMeta}}}()
        for (proj, merged, dm) in valid
            gdir = dirname(proj.dir)
            push!(get!(groups, gdir, []), (proj, merged, dm))
        end

        # collect data across all groups
        run_dfs      = DataFrame[]
        run_scols    = Vector{String}[]
        all_scols    = String[]
        group_labels = String[]       # group name per sample
        run_labels   = String[]       # run name per sample (for shape)
        group_names  = String[]

        for (gdir, members) in groups
            gname = basename(gdir)
            push!(group_names, gname)
            for (proj, merged, dm) in members
                src_csv = merged.tables[_source_key(merged)]
                isfile(src_csv) || continue
                df, scols = _cached_read(cache, src_csv, dm)
                isempty(scols) && continue
                rname = basename(proj.dir)
                push!(run_dfs, df)
                push!(run_scols, scols)
                append!(all_scols, scols)
                append!(group_labels, fill(gname, length(scols)))
                append!(run_labels, fill(rname, length(scols)))
            end
        end

        isempty(run_dfs) && return

        # study NMDS
        mat, _, _ = _build_combined_counts(run_dfs, run_scols, default_db_meta.levels)
        nmds_stress = NaN
        if size(mat, 1) >= 3
            coords, stress = _run_nmds(mat, r_lock)
            if !any(isnan, coords)
                nmds_stress = stress
                _plot() do
                    nmds_plot(coords, all_scols, nmds_pdf;
                              colour_by=group_labels, shape_by=run_labels,
                              colour_label="Group", shape_label="Run",
                              stress=stress, subtitle)
                end
                @info "Written: $nmds_pdf"
                log_written(study_dir, nmds_pdf)
            else
                @warn "Study-level NMDS failed"
                pipeline_log(study_dir, "WARN: Study-level NMDS failed")
            end
        else
            @warn "Too few samples ($(size(mat, 1))) for study-level NMDS"
        end

        # study alpha boxplot by group
        all_alpha    = DataFrame[]
        alpha_labels = String[]
        for (gdir, members) in groups
            gname = basename(gdir)
            combined_alpha = DataFrame(sample=String[], richness=Int[],
                                       shannon=Float64[], simpson=Float64[])
            for (_, merged, dm) in members
                src_csv = merged.tables[_source_key(merged)]
                isfile(src_csv) || continue
                df, scols = _cached_read(cache, src_csv, dm)
                isempty(scols) && continue
                append!(combined_alpha, _compute_alpha(df, scols))
            end
            if nrow(combined_alpha) > 0
                push!(all_alpha, combined_alpha)
                push!(alpha_labels, gname)
            end
        end

        if !isempty(all_alpha)
            _plot() do
                alpha_boxplot(all_alpha, alpha_labels, alpha_pdf; subtitle)
            end
            @info "Written: $alpha_pdf"
            log_written(study_dir, alpha_pdf)
        end

        # PERMANOVA (only when metadata.csv exists)
        metadata = load_metadata(study_dir, study_dir)
        if !isnothing(metadata) && size(mat, 1) >= 3
            sample_col = "sample" in names(metadata) ? "sample" :
                         "Sample" in names(metadata) ? "Sample" : nothing
            if !isnothing(sample_col)
                meta_idx = indexin(all_scols, String.(metadata[!, sample_col]))
                if all(!isnothing, meta_idx)
                    meta_matched = metadata[collect(meta_idx), :]
                    perm_result  = _run_permanova(mat, meta_matched, r_lock)
                    if !isnothing(perm_result)
                        mkpath(dirname(perm_txt))
                        write(perm_txt, perm_result)
                        @info "Written: $perm_txt"
                        log_written(study_dir, perm_txt)
                    end
                else
                    @warn "Not all samples found in metadata.csv - skipping PERMANOVA"
                end
            else
                @warn "metadata.csv has no 'sample' column - skipping PERMANOVA"
            end
        end

        # dual analysis reports
        report_sections = Pair{String, String}[]

        # Group overview.
        buf = IOBuffer()
        println(buf, "Groups: ", join(group_names, ", "))
        println(buf, "Total samples: ", length(all_scols))
        push!(report_sections, "Overview" => String(take!(buf)))

        # NMDS info.
        if !isnan(nmds_stress)
            quality = nmds_stress < 0.2 ? "Good: below 0.2 threshold" : "Poor: above 0.2 threshold"
            push!(report_sections, "NMDS" =>
                  "Stress: $(round(nmds_stress; digits=3)) ($quality)")
        end

        # Alpha diversity by group.
        if !isempty(all_alpha)
            buf = IOBuffer()
            print(buf, rpad("Group", 20))
            println(buf, rpad("Samples", 10), rpad("Mean richness", 16),
                    rpad("Mean Shannon", 16), "Mean Simpson")
            for (adf, gname) in zip(all_alpha, alpha_labels)
                n = nrow(adf)
                mr = round(mean(adf.richness); digits=0)
                ms = round(mean(adf.shannon); digits=3)
                mp = round(mean(adf.simpson); digits=3)
                print(buf, rpad(gname, 20))
                if n > 1
                    sr = round(std(adf.richness); digits=0)
                    ss = round(std(adf.shannon); digits=2)
                    sp = round(std(adf.simpson); digits=2)
                    println(buf, rpad(string(n), 10),
                            rpad("$(Int(mr)) +/- $(Int(sr))", 16),
                            rpad("$ms +/- $ss", 16),
                            "$mp +/- $sp")
                else
                    println(buf, rpad(string(n), 10),
                            rpad(string(Int(mr)), 16),
                            rpad(string(ms), 16),
                            string(mp))
                end
            end
            push!(report_sections, "Alpha Diversity by Group" => String(take!(buf)))
        end

        # Write one report per source per method.
        for method in _TAX_METHODS
            msrc_keys = _method_source_keys(valid[1][2], method)
            for src in msrc_keys
                src_run_dfs  = DataFrame[]
                src_all_scols = String[]
                # Determine the table key to look up
                table_key = src
                for (gdir_inner, members_inner) in groups
                    for (proj, merged, dm) in members_inner
                        csv = get(merged.tables, table_key, "")
                        # For "merged" under dada2, fall back to "merged"
                        if isempty(csv) && method == "dada2" && src == "merged"
                            csv = get(merged.tables, "merged", "")
                        end
                        (isempty(csv) || !isfile(csv)) && continue
                        df, sc = _cached_read(cache, csv, dm)
                        isempty(sc) && continue
                        push!(src_run_dfs, df)
                        append!(src_all_scols, sc)
                    end
                end
                isempty(src_run_dfs) && continue
                method == "dada2" && !_has_dada2(src_run_dfs[1], default_db_meta.levels) && continue

                view_dfs = [_method_df(df, method, default_db_meta.levels) for df in src_run_dfs]
                combined = reduce((a, b) -> vcat(a, b; cols=:union), view_dfs)
                sd = _method_source_dirname(src)
                sl = sd == "unfiltered" ? "unfiltered" : sd
                sections = copy(report_sections)
                for (rank_name, rank_col) in report_ranks
                    section = _top_taxa_section(combined, src_all_scols, rank_name, rank_col; n=20)
                    !isempty(section) && push!(sections, "Top 20 $(rank_name) ($(sl))" => section)
                end

                rpath = joinpath(analysis_dir, "analysis_report_$(method)_$(sd).txt")
                _write_report(rpath,
                              "Study Report\n  Source: $sl ($method)",
                              sections)
                log_written(study_dir, rpath)
            end
        end
    end

    #  Public: analyse_study (entry point from main.jl)
    """
        analyse_study(projects, merged_results)

    Run group-level and study-level analysis after the per-run `@threads`
    loop has completed.

    Groups projects by `dirname(project.dir)` (= group directory).
    Calls `_analyse_group` for each group with ≥2 runs, then
    `_analyse_study_level` when there are ≥2 groups.
    """
    function analyse_study(projects::Vector{ProjectCtx},
                           merged_results::Vector{<:Union{MergedTables, Nothing}},
                           db_metas::Vector{DatabaseMeta};
                           plot_lock::Union{Nothing,ReentrantLock}=nothing)
        isempty(projects) && return

        r_lock    = ReentrantLock()
        study_dir = projects[1].study_dir

        # Pair projects with results and db_metas, filtering out failures.
        valid = Tuple{ProjectCtx, MergedTables, DatabaseMeta}[
            (projects[i], merged_results[i], db_metas[i])
            for i in eachindex(projects)
            if i <= length(merged_results) && !isnothing(merged_results[i])
        ]
        isempty(valid) && return

        # Group by parent directory.
        groups = Dict{String, Vector{Tuple{ProjectCtx, MergedTables, DatabaseMeta}}}()
        for (proj, merged, dm) in valid
            gdir = dirname(proj.dir)
            push!(get!(groups, gdir, []), (proj, merged, dm))
        end

        # Per-group analysis (skipped for single-run groups).
        # Each group uses the db_meta of its first member (all runs in a group share a DB).
        for (gdir, members) in groups
            if length(members) > 1
                group_db_meta = members[1][3]
                group_members = [(p, m) for (p, m, _) in members]
                _analyse_group(gdir, group_members, r_lock, group_db_meta; plot_lock)
            end
        end

        # Study-level analysis (only meaningful with ≥2 genuine groups).
        if length(groups) >= 2
            study_db_meta = valid[1][3]
            _analyse_study_level(study_dir, valid, r_lock, study_db_meta; plot_lock)
        end
    end

end
