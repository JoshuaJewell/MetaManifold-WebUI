# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

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

    # Union source keys across multiple MergedTables (for study/group levels where
    # different runs may carry different filter sets).
    function _all_method_source_keys(mergeds::Vector{MergedTables}, method::String)
        seen = Set{String}()
        out  = String[]
        for m in mergeds
            for k in _method_source_keys(m, method)
                if k ∉ seen
                    push!(seen, k)
                    push!(out, k)
                end
            end
        end
        # Ensure "merged" is always last.
        filter!(!=("merged"), out)
        push!(out, "merged")
        return out
    end

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
    function _source_label(key::String)
        key == "merged"     && return "unfiltered"
        key == "merged_otu" && return "unfiltered (OTU)"
        endswith(key, "_otu") ? replace(key, r"_otu$" => " (OTU)") : key
    end

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

    # Total sequence count across all sample columns. Used to skip empty filtered tables.
    _total_seqs(df::DataFrame, scols::Vector{String}) =
        sum(col -> sum(skipmissing(df[!, col]); init=0), scols; init=0)
