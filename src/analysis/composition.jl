# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

    ## Source table selection
    # Pick the primary analysis source from a MergedTables.
    # Currently always "merged" (unfiltered); filter selection can be added here.
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
