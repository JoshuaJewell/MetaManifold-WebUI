module TaxonomyTableTools

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

using CSV, DataFrames, Logging, YAML
using ..PipelineTypes
using ..PipelineGraph
using ..PipelineLog
using ..Config

export merge_taxonomy_counts, filter_table, filter_table_dada2, merge_taxa

    # Import vsearch taxonomy
    function import_vsearch(file::AbstractString)
        rows = Vector{Vector{String}}()
        open(file, "r") do io
            for line in eachline(io)
                clean = replace(chomp(line), "\r" => "")
                isempty(clean) && continue
                push!(rows, split(clean, '\t'))
            end
        end
        return rows
    end

    function vsearch_to_df(imported_list::Vector{Vector{String}})
        qseqid = String[]
        pident = Float64[]
        sseqtax = String[]
        for rec in imported_list
            length(rec) < 3 && continue
            push!(qseqid, rec[1])
            push!(sseqtax, rec[2])
            push!(pident, parse(Float64, rec[3]))
        end
        return DataFrame(SeqName=qseqid, Pident=pident, sseqtax=sseqtax)
    end

    function build_taxonomy_table_rows(df::DataFrame, db_meta::DatabaseMeta)
        levels = db_meta.levels
        if db_meta.vsearch_format == "pr2"
            # PR2: pipe-separated with metadata prefix
            meta_cols = ["SeqName", "Pident", "Accession", "rRNA", "Organellum", "specimen"]
            header = vcat(meta_cols, levels)
            n_meta = length(meta_cols) - 2  # fields parsed from sseqtax (excl SeqName, Pident)
        else
            # Generic: semicolon-separated, mapped positionally to levels
            header = vcat(["SeqName", "Pident"], levels)
            n_meta = 0
        end

        out_rows = Vector{Vector{String}}(undef, nrow(df))

        for (i, r) in enumerate(eachrow(df))
            data = fill("", length(header))
            data[1] = string(r.SeqName)
            data[2] = string(r.Pident)

            tax_str = r.sseqtax
            if db_meta.vsearch_format == "pr2" && occursin('|', tax_str)
                # Full PR2 format: Accession|rRNA|Organellum|specimen|Domain|...|Species
                parts = split(tax_str, '|')
                for j in 1:min(length(parts), length(header) - 2)
                    data[j + 2] = parts[j]
                end
            elseif db_meta.vsearch_format == "pr2"
                # PR2 taxonomy-only: Domain;Supergroup;...;Species[;]
                parts = filter(!isempty, split(tax_str, ';'))
                tax_start = length(header) - length(levels) + 1
                for j in 1:min(length(parts), length(levels))
                    data[tax_start + j - 1] = parts[j]
                end
            else
                # Generic: semicolon-separated, mapped positionally to levels
                parts = filter(!isempty, split(tax_str, ';'))
                tax_start = 3
                for j in 1:min(length(parts), length(levels))
                    data[tax_start + j - 1] = String(strip(parts[j]))
                end
            end

            out_rows[i] = data
        end
        return header, out_rows
    end

    # Sorting helper for DADA2 seq IDs (seq1, seq2, ...).
    function seqnum(x)
        if ismissing(x)
            return typemax(Int)
        end
        m = match(r"seq(\d+)", String(x))
        return m === nothing ? typemax(Int) : parse(Int, m.captures[1])
    end

    """
        merge_taxonomy_counts(taxonomy_vsearch_path, counts_csv_path, db_meta)

    Reads taxonomy TSV from vsearch and counts CSV from DADA2 and returns a merged DataFrame.

    ## Arguments:
    - `taxonomy_vsearch_path`: Path to vsearch output file.
    - `counts_csv_path`: Path to DADA2 idtax output file.
    - `db_meta`: DatabaseMeta with levels, corrections, and vsearch_format.
    """
    function merge_taxonomy_counts(
        taxonomy_vsearch_path,
        counts_csv_path,
        db_meta::DatabaseMeta;
        bootstraps_path=nothing)
        
        @info("Merging taxonomy counts.")

        # build taxonomy dataframe - FIXED VERSION
        imported = import_vsearch(taxonomy_vsearch_path)
        df_tax_raw = vsearch_to_df(imported)
        header, rows = build_taxonomy_table_rows(df_tax_raw, db_meta)
        
        # Create DataFrame properly by specifying all columns first
        df_taxonomy = DataFrame([Symbol(h) => String[] for h in header])
        
        # Append each row properly
        for r in rows
            push!(df_taxonomy, r)
        end

        # read counts CSV
        df_counts = CSV.read(counts_csv_path, DataFrame)

        # detect columns: seq id and optional sequence column
        seq_id_col = nothing
        sequence_col = nothing

        for col in names(df_counts)
            sample = collect(skipmissing(df_counts[!, col]))
            isempty(sample) && continue
            first_val = string(first(sample))

            if occursin(r"^seq\d+$", first_val)
                seq_id_col = col
            end

            if length(first_val) > 50 && occursin(r"^[ACGTN]+$", first_val)
                sequence_col = col
            end
        end

        # Fallback: if no data rows, detect seq-ID column by name
        if isnothing(seq_id_col)
            "SeqName" in names(df_counts) && (seq_id_col = "SeqName")
        end
        if isnothing(seq_id_col)
            error("Cannot find a column with seq IDs (seq1, seq2, ...) in counts file!")
        end

        df_counts_prepared = copy(df_counts)

        # If sequence column exists and isn't named "sequence", rename or drop accordingly
        if !isnothing(sequence_col) && sequence_col != "sequence"
            if "sequence" in names(df_counts_prepared)
                select!(df_counts_prepared, Not(sequence_col))
            else
                rename!(df_counts_prepared, sequence_col => "sequence")
            end
        end

        # Ensure SeqName column exists and refers to seq id column
        if seq_id_col != "SeqName"
            if "SeqName" in names(df_counts_prepared)
                if isnothing(sequence_col)
                    rename!(df_counts_prepared, "SeqName" => "sequence")
                else
                    select!(df_counts_prepared, Not("SeqName"))
                end
            end
            rename!(df_counts_prepared, seq_id_col => "SeqName")
        end

        # Rename any columns in the DADA2 file that clash with vsearch taxonomy
        # columns (e.g. Domain, Supergroup, ... when counts_csv_path is taxonomy.csv)
        overlap = filter(!=(Symbol("SeqName")),
                         intersect(Symbol.(names(df_taxonomy)),
                                   Symbol.(names(df_counts_prepared))))
        if !isempty(overlap)
            rename!(df_counts_prepared,
                    [String(c) => String(c) * "_dada2" for c in overlap])
        end

        # Full outer join to retain ASVs found by either method
        merged_df = outerjoin(df_taxonomy, df_counts_prepared, on="SeqName")
        sort!(merged_df, "SeqName", by=seqnum)

        # Join bootstrap confidence values if a bootstraps file was provided.
        if !isnothing(bootstraps_path) && isfile(bootstraps_path) && filesize(bootstraps_path) > 0
            df_boot = CSV.read(bootstraps_path, DataFrame)
            # Keep only SeqName and *_boot columns; drop Sequence (already in merged_df).
            boot_cols = filter(c -> c == "SeqName" || endswith(c, "_boot"), names(df_boot))
            select!(df_boot, boot_cols)
            merged_df = leftjoin(merged_df, df_boot, on="SeqName")
            @info "Joined $(length(boot_cols)-1) bootstrap columns from $bootstraps_path"
        end

        # Fill NA/missing values at lower ranks with parent_X convention.
        # e.g. if Family=Muribaculaceae and Genus=NA, Genus becomes "Muribaculaceae_X".
        # Applied to both plain rank columns and _dada2 suffixed columns.
        for suffix in ("", "_dada2")
            cols = [l * suffix for l in db_meta.levels if (l * suffix) in names(merged_df)]
            for i in 2:length(cols)
                col = cols[i]
                parent_col = cols[i-1]
                merged_df[!, col] = [
                    begin
                        v = merged_df[rn, col]
                        is_na = ismissing(v) || strip(string(v)) == "NA" || strip(string(v)) == ""
                        if is_na
                            p = merged_df[rn, parent_col]
                            (ismissing(p) || strip(string(p)) == "NA" || strip(string(p)) == "") ?
                                v : string(strip(string(p))) * "_X"
                        else
                            v
                        end
                    end
                    for rn in 1:nrow(merged_df)
                ]
            end
        end

        # Apply database-specific taxonomy corrections from config.
        for corr in db_meta.corrections
            src_col = get(corr, "source", nothing)
            tgt_col = get(corr, "target", nothing)
            vals    = get(corr, "values", Dict())
            (isnothing(src_col) || isnothing(tgt_col)) && continue
            src_col, tgt_col = string(src_col), string(tgt_col)
            if src_col in names(merged_df) && tgt_col in names(merged_df) && !isempty(vals)
                mapping = Dict(string(k) => string(v) for (k,v) in vals)
                merged_df[!, tgt_col] = [
                    ismissing(d) ? s : get(mapping, string(d), ismissing(s) ? s : string(s))
                    for (d, s) in zip(merged_df[!, src_col], merged_df[!, tgt_col])
                ]
            end
        end

        @info "Merge complete. $(nrow(merged_df)) rows."

        return merged_df
    end

    """
        filter_table(merged_df, filters_yaml_path)

    Apply taxonomy filtering to a merged DataFrame using rules from a YAML file.

    ## Arguments:
    - `merged_df`: Input DataFrame with taxonomy columns.
    - `filters_yaml_path`: Path to YAML file defining mappings, filters, and remove_empty rules.

    ## YAML structure:
    ``` YAML
        mappings:
          - source_column: Division
            target_column: Supergroup
            values: { Rhizaria: Rhizaria, ... }
        filters:
          - column: Domain
            pattern: Bacteria
            action: exclude    # exclude (default) | keep
            regex: false       # false (default) = substring match
        remove_empty:
          - Domain
    ```

    A flat `mappings` dict and `remove_empty_domain` bool are auto-detected as legacy format.
    """
    function _apply_single_mapping!(df::DataFrame, src_col::String, tgt_col::String, mapping)
        isempty(mapping) && return
        m = Dict(string(k) => string(v) for (k,v) in mapping)
        if src_col in names(df) && tgt_col in names(df)
            df[!, tgt_col] = [
                ismissing(s) ? t : get(m, string(s), t)
                for (s, t) in zip(df[!, src_col], df[!, tgt_col])
            ]
        elseif !isempty(m)
            missing_cols = filter(c -> c ∉ names(df), [src_col, tgt_col])
            @warn "Column remapping skipped: $(join(missing_cols, ", ")) not in data. " *
                  "Available: $(join(names(df), ", "))"
        end
    end

    function _apply_mappings!(df::DataFrame, mapping_config)
        if mapping_config isa AbstractDict
            # Legacy: flat dict = Division -> Supergroup
            _apply_single_mapping!(df, "Division", "Supergroup", mapping_config)
        elseif mapping_config isa Vector
            for entry in mapping_config
                src = get(entry, "source_column", nothing)
                tgt = get(entry, "target_column", nothing)
                vals = get(entry, "values", Dict())
                (isnothing(src) || isnothing(tgt)) && continue
                _apply_single_mapping!(df, string(src), string(tgt), vals)
            end
        end
    end

    function filter_table(
        merged_df::DataFrame,
        filters_yaml_path::String
    )
        config = YAML.load_file(filters_yaml_path)
        df = deepcopy(merged_df)

        @info "Filtering table using configuration from $filters_yaml_path"

        # mappings (auto-detect old vs new format)
        _apply_mappings!(df, get(config, "mappings", Dict()))

        # remove_empty (backwards compat with remove_empty_domain)
        remove_cols = if haskey(config, "remove_empty")
            val = config["remove_empty"]
            val isa Vector ? string.(val) : [string(val)]
        elseif get(config, "remove_empty_domain", false)
            ["Domain"]
        else
            String[]
        end

        for col in remove_cols
            if col in names(df)
                df = filter(row ->
                    !ismissing(row[col]) &&
                    !isempty(strip(string(row[col]))) &&
                    lowercase(strip(string(row[col]))) != "blank",
                    df)
            else
                @warn "remove_empty: column '$col' not in data. Available: $(join(names(df), ", "))"
            end
        end

        # filters (with action + regex support)
        raw_filters = get(config, "filters", [])
        @assert raw_filters isa Vector "filters must be a list in YAML"

        for item in raw_filters
            @assert item isa Dict "Each filter must be a key-value map"
            colname   = string(get(item, "column", ""))
            pattern   = string(get(item, "pattern", ""))
            action    = lowercase(string(get(item, "action", "exclude")))
            use_regex = get(item, "regex", false) == true

            if !(colname in names(df))
                @warn "Skipping filter: column '$colname' not present in data. " *
                      "Available: $(join(names(df), ", "))"
                continue
            end

            matcher = use_regex ? (val -> occursin(Regex(pattern), string(val))) :
                                  (val -> occursin(pattern, string(val)))

            if action == "keep"
                df = filter(row -> ismissing(row[colname]) || matcher(row[colname]), df)
            else  # exclude (default)
                df = filter(row -> ismissing(row[colname]) || !matcher(row[colname]), df)
            end
        end

        @info "Filtering complete. $(nrow(df))/$(nrow(merged_df)) rows retained."

        return df
    end

    """
        _with_dada2_as_primary(df, levels) -> DataFrame

    Return a copy of `df` where `_dada2` taxonomy columns are renamed to
    standard names (and original vsearch columns get a `_vsearch` suffix).
    This lets `filter_table` operate on DADA2 taxonomy using unchanged logic.
    """
    function _with_dada2_as_primary(df::DataFrame, levels::Vector{String})
        out = copy(df)
        for col in levels
            dada2_col = col * "_dada2"
            dada2_col in names(out) || continue
            if col in names(out)
                rename!(out, col => col * "_vsearch")
            end
            rename!(out, dada2_col => col)
        end
        # Also swap Pident if a DADA2 equivalent exists
        if "Pident_dada2" in names(out) && "Pident" in names(out)
            rename!(out, "Pident" => "Pident_vsearch", "Pident_dada2" => "Pident")
        end
        return out
    end

    """
        _restore_from_dada2_primary(df, levels) -> DataFrame

    Reverse of `_with_dada2_as_primary`: rename standard columns back to
    `_dada2` and `_vsearch` columns back to standard names.
    """
    function _restore_from_dada2_primary(df::DataFrame, levels::Vector{String})
        out = copy(df)
        for col in levels
            vsearch_col = col * "_vsearch"
            # Current standard name holds DADA2 values - rename to _dada2
            if col in names(out)
                rename!(out, col => col * "_dada2")
            end
            # Restore vsearch column
            if vsearch_col in names(out)
                rename!(out, vsearch_col => col)
            end
        end
        if "Pident_vsearch" in names(out)
            if "Pident" in names(out)
                rename!(out, "Pident" => "Pident_dada2")
            end
            rename!(out, "Pident_vsearch" => "Pident")
        end
        return out
    end

    """
        filter_table_dada2(merged_df, filters_yaml_path, levels)

    Apply taxonomy filtering using DADA2 columns instead of vsearch columns.
    Temporarily swaps column names so `filter_table` logic applies to `_dada2` columns.
    """
    function filter_table_dada2(merged_df::DataFrame, filters_yaml_path::String,
                                levels::Vector{String})
        swapped = _with_dada2_as_primary(merged_df, levels)
        filtered = filter_table(swapped, filters_yaml_path)
        return _restore_from_dada2_primary(filtered, levels)
    end

    function merge_taxa(project::ProjectCtx, source::ASVResult, tax::TaxonomyHits,
                        db_meta::DatabaseMeta)
        config_path = write_run_config(project)
        cfg         = get(YAML.load_file(config_path), stage_sections(:merge_taxa), Dict())
        filter_list = get(cfg, "filters", [nothing])
        merge_dir   = joinpath(project.dir, "merged")
        hash_file   = joinpath(merge_dir, "config.hash")

        tax_counts_path = joinpath(dirname(source.taxonomy), "tax_counts.csv")
        tax_prefix_     = splitext(basename(source.taxonomy))[1]
        boot_path_      = joinpath(dirname(source.taxonomy), tax_prefix_ * "_bootstraps.csv")
        data_mtime     = max(mtime(tax.tsv), mtime(source.taxonomy),
                             isfile(tax_counts_path) ? mtime(tax_counts_path) : 0.0,
                             isfile(boot_path_)      ? mtime(boot_path_)      : 0.0)
        config_changed = _section_stale(config_path, stage_sections(:merge_taxa), hash_file)
        merged_csv     = joinpath(merge_dir, "merged.csv")

        # Determine what needs to be (re-)computed.
        need_base     = config_changed || !isfile(merged_csv) || mtime(merged_csv) <= data_mtime
        stale_filters = Pair{String,String}[]  # stem => filter_path
        tables        = Dict{String,String}("merged" => merged_csv)

        filter_colours = Dict{String,String}()
        # Pre-compute applicable filter stems (respecting per-filter database restrictions).
        filter_stems = String[]
        for e in filter_list
            isnothing(e) && continue
            fp = joinpath(project.config_dir, "filters", string(e))
            if isfile(fp)
                allowed = get(YAML.load_file(fp), "databases", nothing)
                isnothing(allowed) || db_meta.name ∈ string.(allowed) || continue
            end
            push!(filter_stems, splitext(string(e))[1])
        end

        for entry in filter_list
            isnothing(entry) && continue
            filter_name  = string(entry)
            filter_path  = joinpath(project.config_dir, "filters", filter_name)
            stem         = splitext(filter_name)[1]
            output_csv   = joinpath(merge_dir, stem * ".csv")

            # Skip filters that declare a databases list not including the active DB.
            if isfile(filter_path)
                fcfg = YAML.load_file(filter_path)
                allowed_dbs = get(fcfg, "databases", nothing)
                if !isnothing(allowed_dbs) && db_meta.name ∉ string.(allowed_dbs)
                    @info "Skipping merge_taxa filter '$stem': not applicable to database '$(db_meta.name)'"
                    continue
                end
                if haskey(fcfg, "colour")
                    filter_colours[stem] = string(fcfg["colour"])
                end
            end

            tables[stem] = output_csv
            # Pre-populate the DADA2-filtered path so it is always present in the
            # returned MergedTables, even when the function returns early via the
            # skip guard below (the CSV may already exist from a previous run).
            tables[stem * "_dada2"] = joinpath(merge_dir, stem * "_dada2" * ".csv")

            if config_changed || !isfile(output_csv) ||
               mtime(output_csv) <= max(data_mtime, mtime(filter_path))
                push!(stale_filters, stem => filter_path)
            else
                @info "Skipping merge_taxa filter '$stem': $output_csv up to date"
            end
        end

        if !need_base && isempty(stale_filters)
            @info "Skipping merge_taxa: all outputs up to date in $merge_dir"
            return MergedTables(tables, filter_stems, filter_colours)
        end

        mkpath(merge_dir)
        counts_path = (isfile(tax_counts_path) && filesize(tax_counts_path) > 0) ?
                      tax_counts_path : source.taxonomy
        tables_dir = dirname(source.taxonomy)
        tax_prefix = splitext(basename(source.taxonomy))[1]   # e.g. "taxonomy"
        boot_path  = joinpath(tables_dir, tax_prefix * "_bootstraps.csv")
        df = need_base ? merge_taxonomy_counts(tax.tsv, counts_path, db_meta;
                             bootstraps_path=boot_path) :
                         CSV.read(merged_csv, DataFrame)

        if need_base
            CSV.write(merged_csv, df)
            @info "Written: $merged_csv"
            pipeline_log(project, "Merge complete. $(nrow(df)) ASVs.")
            log_written(project, merged_csv)
        end

        for (stem, filter_path) in stale_filters
            output_csv = tables[stem]
            filtered = filter_table(df, filter_path)
            CSV.write(output_csv, filtered)
            @info "Written: $output_csv"
            pipeline_log(project, "Filter '$stem': $(nrow(filtered))/$(nrow(df)) rows retained.")
            log_written(project, output_csv)
        end

        # DADA2 filtered CSVs (same filters, applied to _dada2 columns)
        # Check whether any _dada2 taxonomy columns exist before attempting.
        has_dada2 = any(endswith(c, "_dada2") for c in names(df)
                        if any(startswith(c, t) for t in db_meta.levels))

        if has_dada2
            for entry in filter_list
                isnothing(entry) && continue
                filter_name  = string(entry)
                filter_path  = joinpath(project.config_dir, "filters", filter_name)
                stem         = splitext(filter_name)[1]
                dada2_stem   = stem * "_dada2"
                output_csv   = joinpath(merge_dir, dada2_stem * ".csv")

                # Skip filters that declare a databases list not including the active DB.
                if isfile(filter_path)
                    allowed_dbs = get(YAML.load_file(filter_path), "databases", nothing)
                    if !isnothing(allowed_dbs) && db_meta.name ∉ string.(allowed_dbs)
                        @info "Skipping merge_taxa DADA2 filter '$stem': not applicable to database '$(db_meta.name)'"
                        continue
                    end
                end

                if !config_changed && isfile(output_csv) &&
                   mtime(output_csv) > max(data_mtime, mtime(filter_path))
                    @info "Skipping merge_taxa DADA2 filter '$stem': $output_csv up to date"
                    continue
                end

                filtered = filter_table_dada2(df, filter_path, db_meta.levels)
                CSV.write(output_csv, filtered)
                @info "Written: $output_csv (DADA2 filter)"
                pipeline_log(project, "Filter '$stem' (DADA2): $(nrow(filtered))/$(nrow(df)) rows retained.")
                log_written(project, output_csv)
            end
        end

        _write_section_hash(config_path, stage_sections(:merge_taxa), hash_file)
        return MergedTables(tables, filter_stems, filter_colours)
    end

end