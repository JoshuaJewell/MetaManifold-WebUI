module PipelinePlotsPlotly

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Plotly JSON figure generation for the WebUI.
#
# Each public function takes the same logical data as its PipelinePlots
# counterpart and writes a Plotly-compatible JSON figure spec to the given
# path. No CairoMakie dependency - pure Julia + JSON3.
#
# Plotly figure JSON structure:
#   { "data": [ <traces> ], "layout": { <layout> } }

export taxa_bar_chart, filter_composition_plot, alpha_diversity_plot,
       pipeline_stats_plot, nmds_plot, alpha_boxplot, group_comparison_chart

    using DataFrames, JSON3

    ## Color palette

    # Wong palette (colour-blindness friendly) as hex, extended via golden-angle
    # hue rotation - same logical palette as PipelinePlots._palette.
    function _palette_hex(n::Int)::Vector{String}
        base = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7"]
        n <= length(base) && return base[1:n]
        colors = copy(base)
        for i in (length(base) + 1):n
            hue = mod(i * 137.508, 360.0)
            s, v = 0.65, 0.85
            c = v * s
            x = c * (1.0 - abs(mod(hue / 60.0, 2.0) - 1.0))
            m = v - c
            r, g, b = hue < 60  ? (c + m, x + m, m) :
                    hue < 120 ? (x + m, c + m, m) :
                    hue < 180 ? (m, c + m, x + m) :
                    hue < 240 ? (m, x + m, c + m) :
                    hue < 300 ? (x + m, m, c + m) :
                                (c + m, m, x + m)
            push!(colors, "#" * lpad(string(round(Int, r * 255), base = 16), 2, '0') *
                                lpad(string(round(Int, g * 255), base = 16), 2, '0') *
                                lpad(string(round(Int, b * 255), base = 16), 2, '0'))
        end
        colors[1:n]
    end

    const _GREY_HEX = "#B3B3B3"

    function _apply_grey!(colors::Vector{String}, labels::Vector{String})
        for (i, l) in enumerate(labels)
            l in ("Unclassified", "Other") && (colors[i] = _GREY_HEX)
        end
    end

    ## Display-name helper (same logic as PipelinePlots._display_name)
    const _FASTQ_SUFFIX = r"_R[12](?:_filt|_trimmed)?\.fastq\.gz$"
    _display_name(s::String)  = replace(s, _FASTQ_SUFFIX => "")
    _display_names(v::Vector{String}) = String[_display_name(s) for s in v]

    ## Shared data helpers (same logic as PipelinePlots, duplicated for module independence)

    function _lowest_rank_label(row; rank_order::Vector{String})
        label = "Unclassified"
        for rank in rank_order
            if hasproperty(row, Symbol(rank))
                val = row[Symbol(rank)]
                if !ismissing(val) && !isempty(strip(string(val)))
                    label = string(val)
                end
            end
        end
        label
    end

    function _rank_label(row, rank::String)
        sym = Symbol(rank)
        if hasproperty(row, sym)
            val = row[sym]
            if !ismissing(val) && !isempty(strip(string(val)))
                return string(val)
            end
        end
        "Unclassified"
    end

    function _label_rows(df::DataFrame; rank = nothing,
                        rank_order::Vector{String})::Vector{String}
        if isnothing(rank)
            String[_lowest_rank_label(row; rank_order) for row in eachrow(df)]
        else
            String[_rank_label(row, rank) for row in eachrow(df)]
        end
    end

    function _aggregate_taxa(df::DataFrame, sample_cols::Vector{String},
                            labels::AbstractVector{<:AbstractString};
                            top_n::Int = 15, relative::Bool = true)
        unique_labels = String.(unique(labels))
        label_idx = Dict(l => i for (i, l) in enumerate(unique_labels))
        n_labels  = length(unique_labels)
        n_samples = length(sample_cols)

        counts = zeros(Float64, n_labels, n_samples)
        for (ri, row) in enumerate(eachrow(df))
            li = label_idx[labels[ri]]
            for (j, col) in enumerate(sample_cols)
                v = row[Symbol(col)]
                counts[li, j] += ismissing(v) ? 0.0 :
                    Float64(v isa AbstractString ? parse(Float64, v) : v)
            end
        end

        totals = vec(sum(counts, dims = 2))
        order  = sortperm(totals, rev = true)

        if n_labels > top_n
            keep   = order[1:top_n]
            other  = vec(sum(counts[order[(top_n + 1):end], :], dims = 1))
            counts = vcat(counts[keep, :], other')
            final_labels = vcat(unique_labels[keep], ["Other"])
        else
            counts = counts[order, :]
            final_labels = unique_labels[order]
        end

        if relative
            col_sums = vec(sum(counts, dims = 1))
            for j in 1:n_samples
                col_sums[j] > 0 && (counts[:, j] ./= col_sums[j])
            end
        end

        counts, final_labels
    end

    const _STAGE_LABELS = Dict(
        "input"     => "Input",
        "filtered"  => "Filtered",
        "denoisedF" => "Denoised (F)",
        "denoisedR" => "Denoised (R)",
        "merged"    => "Merged",
        "nochim"    => "Chimera-free",
    )
    const _STAGE_ORDER = ["input", "filtered", "denoisedF", "denoisedR", "merged", "nochim"]

    function _stage_label(col::String)
        haskey(_STAGE_LABELS, col) && return _STAGE_LABELS[col]
        startswith(col, "reads_") && return "After " * replace(col[7:end], "_" => " ")
        col
    end

    ## Title helper: appends subtitle as a small superscript if provided.
    _title(main::String, subtitle) =
        isnothing(subtitle) ? main : "$main<br><sup>$subtitle</sup>"

    ## File writer
    _write_json(path::String, fig::Dict) = write(path, JSON3.write(fig))

    ## Public: taxa_bar_chart
    """
        taxa_bar_chart(df, sample_cols, output_path; top_n, rank, relative,
                    subtitle, rank_order)

    Stacked bar chart of taxonomic composition per sample, written as Plotly JSON.
    """
    function taxa_bar_chart(df::DataFrame, sample_cols::Vector{String},
                            output_path::String;
                            top_n::Int = 15, rank = nothing,
                            relative::Bool = true,
                            subtitle = nothing,
                            rank_order::Vector{String})
        (nrow(df) == 0 || isempty(sample_cols)) && return nothing

        labels = _label_rows(df; rank, rank_order)
        counts, final_labels = _aggregate_taxa(df, sample_cols, labels; top_n, relative)
        n_taxa    = length(final_labels)
        colors    = _palette_hex(n_taxa)
        _apply_grey!(colors, final_labels)
        display_names = _display_names(sample_cols)
        ylabel = relative ? "Relative abundance" : "Read count"
        title  = relative ? "Taxonomic composition" : "Taxonomic composition (absolute)"

        traces = [Dict{String,Any}(
            "type"   => "bar",
            "name"   => final_labels[i],
            "x"      => display_names,
            "y"      => collect(counts[i, :]),
            "marker" => Dict("color" => colors[i]),
        ) for i in 1:n_taxa]

        layout = Dict{String,Any}(
            "barmode" => "stack",
            "title"   => Dict("text" => _title(title, subtitle)),
            "xaxis"   => Dict("title" => "Sample", "tickangle" => -45),
            "yaxis"   => Dict{String,Any}("title" => ylabel),
            "legend"  => Dict("traceorder" => "normal"),
        )
        relative && (layout["yaxis"]["range"] = [0, 1])

        mkpath(dirname(output_path))
        _write_json(output_path, Dict("data" => traces, "layout" => layout))
        return nothing
    end

    ## Public: filter_composition_plot
    """
        filter_composition_plot(filter_totals, sample_names, output_path;
                                subtitle, colour_overrides)

    Stacked bar chart of per-sample read counts split by biological-group filter.
    """
    function filter_composition_plot(filter_totals::Dict{String, Vector{Float64}},
                                    sample_names::Vector{String},
                                    output_path::String;
                                    subtitle = nothing,
                                    colour_overrides::Dict{String,String} = Dict{String,String}())
        isempty(filter_totals) && return nothing

        display_names = _display_names(sample_names)
        filter_names  = sort(collect(keys(filter_totals)))
        if "Unclassified" in filter_names
            filter!(!=("Unclassified"), filter_names)
            push!(filter_names, "Unclassified")
        end
        n_filters = length(filter_names)
        colors    = _palette_hex(n_filters)
        _apply_grey!(colors, filter_names)

        for (i, fname) in enumerate(filter_names)
            haskey(colour_overrides, fname) && (colors[i] = colour_overrides[fname])
        end

        traces = [Dict{String,Any}(
            "type"   => "bar",
            "name"   => filter_names[i],
            "x"      => display_names,
            "y"      => filter_totals[filter_names[i]],
            "marker" => Dict("color" => colors[i]),
        ) for i in 1:n_filters]

        layout = Dict{String,Any}(
            "barmode" => "stack",
            "title"   => Dict("text" => _title("Filter composition", subtitle)),
            "xaxis"   => Dict("title" => "Sample", "tickangle" => -45),
            "yaxis"   => Dict("title" => "Read count"),
        )

        mkpath(dirname(output_path))
        _write_json(output_path, Dict("data" => traces, "layout" => layout))
        return nothing
    end

    ## Public: alpha_diversity_plot
    """
        alpha_diversity_plot(alpha_df, output_path; subtitle)

    Three-panel bar chart (richness, Shannon, Simpson) as Plotly JSON.
    `alpha_df` must have columns :sample, :richness, :shannon, :simpson.
    """
    function alpha_diversity_plot(alpha_df::DataFrame, output_path::String;
                                subtitle = nothing)
        nrow(alpha_df) == 0 && return nothing
        samples = _display_names(String.(alpha_df.sample))
        color   = _palette_hex(1)[1]

        metrics = [
            ("y",  "y",  "Richness (observed ASVs)", Float64.(alpha_df.richness)),
            ("y2", "x2", "Shannon index",            Float64.(alpha_df.shannon)),
            ("y3", "x3", "Simpson index",            Float64.(alpha_df.simpson)),
        ]

        traces = [Dict{String,Any}(
            "type"       => "bar",
            "x"          => samples,
            "y"          => vals,
            "xaxis"      => xax,
            "yaxis"      => yax,
            "showlegend" => false,
            "marker"     => Dict("color" => color),
        ) for (yax, xax, _, vals) in metrics]

        layout = Dict{String,Any}(
            "title"  => Dict("text" => _title("Alpha diversity", subtitle)),
            "grid"   => Dict("rows" => 3, "columns" => 1, "pattern" => "independent"),
            "yaxis"  => Dict("title" => "Richness (observed ASVs)"),
            "yaxis2" => Dict("title" => "Shannon index"),
            "yaxis3" => Dict("title" => "Simpson index"),
            "xaxis3" => Dict("title" => "Sample", "tickangle" => -45),
            "xaxis"  => Dict("tickangle" => -45),
            "xaxis2" => Dict("tickangle" => -45),
        )

        mkpath(dirname(output_path))
        _write_json(output_path, Dict("data" => traces, "layout" => layout))
        return nothing
    end

    ## Public: pipeline_stats_plot
    """
        pipeline_stats_plot(stats_df, output_path; subtitle)

    Grouped bar chart: read counts at each pipeline stage, one bar group per sample.
    """
    function pipeline_stats_plot(stats_df::DataFrame, output_path::String;
                                subtitle = nothing)
        nrow(stats_df) == 0 && return nothing

        all_cols   = names(stats_df)
        stage_cols = [c for c in _STAGE_ORDER if c in all_cols]
        reads_cols = sort([c for c in all_cols if startswith(c, "reads_")])
        append!(stage_cols, reads_cols)
        isempty(stage_cols) && return nothing

        sample_names  = _display_names(String.(stats_df.sample))
        stage_labels  = [_stage_label(c) for c in stage_cols]
        colors        = _palette_hex(length(sample_names))

        traces = [Dict{String,Any}(
            "type"   => "bar",
            "name"   => sample_names[si],
            "x"      => stage_labels,
            "y"      => [let v = stats_df[si, Symbol(c)]; ismissing(v) ? 0.0 : Float64(v) end
                        for c in stage_cols],
            "marker" => Dict("color" => colors[si]),
        ) for si in eachindex(sample_names)]

        layout = Dict{String,Any}(
            "barmode" => "group",
            "title"   => Dict("text" => _title("Reads through pipeline stages", subtitle)),
            "xaxis"   => Dict("title" => "Pipeline stage", "tickangle" => -45),
            "yaxis"   => Dict("title" => "Read count"),
        )

        mkpath(dirname(output_path))
        _write_json(output_path, Dict("data" => traces, "layout" => layout))
        return nothing
    end

    ## Public: nmds_plot
    """
        nmds_plot(coords, labels, output_path; colour_by, shape_by,
                colour_label, shape_label, stress, subtitle)

    NMDS ordination scatter plot. `coords` is an nx2 Float64 matrix.
    """
    function nmds_plot(coords::Matrix{Float64}, labels::Vector{String},
                    output_path::String;
                    colour_by = nothing,
                    shape_by  = nothing,
                    colour_label::String = "Group",
                    shape_label::String  = "Group",
                    stress = nothing,
                    subtitle = nothing)
        n = size(coords, 1)
        n == 0 && return nothing

        cgroups = isnothing(colour_by) ? ["all"] : unique(colour_by)
        colors  = _palette_hex(length(cgroups))
        cb      = isnothing(colour_by) ? fill("all", n) : colour_by

        traces = map(enumerate(cgroups)) do (ci, cg)
            mask = cb .== cg
            Dict{String,Any}(
                "type"   => "scatter",
                "mode"   => "markers",
                "name"   => cg,
                "x"      => coords[mask, 1],
                "y"      => coords[mask, 2],
                "text"   => labels[mask],
                "marker" => Dict("color" => colors[ci], "size" => 10),
            )
        end

        annotations = Any[]
        if !isnothing(stress)
            push!(annotations, Dict{String,Any}(
                "text"    => "stress = $(round(stress; digits = 3))",
                "xref"    => "paper", "yref" => "paper",
                "x"       => 0.02,    "y"    => 0.98,
                "xanchor" => "left",  "yanchor" => "top",
                "showarrow" => false,
                "font"    => Dict("size" => 12),
            ))
        end

        layout = Dict{String,Any}(
            "title"       => Dict("text" => _title("NMDS ordination", subtitle)),
            "xaxis"       => Dict("title" => "NMDS1"),
            "yaxis"       => Dict("title" => "NMDS2"),
            "annotations" => annotations,
        )

        mkpath(dirname(output_path))
        _write_json(output_path, Dict("data" => traces, "layout" => layout))
        return nothing
    end

    ## Public: alpha_boxplot
    """
        alpha_boxplot(all_alpha, group_labels, output_path; subtitle)

    Three-panel boxplot (richness, Shannon, Simpson) comparing groups.
    """
    function alpha_boxplot(all_alpha::Vector{DataFrame},
                        group_labels::Vector{String},
                        output_path::String;
                        subtitle = nothing)
        isempty(all_alpha) && return nothing
        n_groups = length(group_labels)
        colors   = _palette_hex(n_groups)

        panels = [
            ("y",  "Richness (observed ASVs)", :richness),
            ("y2", "Shannon index",            :shannon),
            ("y3", "Simpson index",            :simpson),
        ]

        traces = Any[]
        for (pi, (yax, _, col)) in enumerate(panels)
            for (gi, (adf, gname)) in enumerate(zip(all_alpha, group_labels))
                push!(traces, Dict{String,Any}(
                    "type"       => "box",
                    "name"       => gname,
                    "y"          => Float64.(adf[!, col]),
                    "yaxis"      => yax,
                    "xaxis"      => pi == 1 ? "x" : "x$(pi)",
                    "marker"     => Dict("color" => colors[gi]),
                    "showlegend" => pi == 1,
                    "legendgroup" => gname,
                ))
            end
        end

        layout = Dict{String,Any}(
            "title"  => Dict("text" => _title("Alpha diversity comparison", subtitle)),
            "grid"   => Dict("rows" => 3, "columns" => 1, "pattern" => "independent"),
            "yaxis"  => Dict("title" => "Richness (observed ASVs)"),
            "yaxis2" => Dict("title" => "Shannon index"),
            "yaxis3" => Dict("title" => "Simpson index"),
        )

        mkpath(dirname(output_path))
        _write_json(output_path, Dict("data" => traces, "layout" => layout))
        return nothing
    end

    ## Public: group_comparison_chart
    """
        group_comparison_chart(dfs_per_group, scols_per_group, group_names,
                            output_path; top_n, rank, relative, subtitle, rank_order)

    Stacked bar chart with one bar per group (reads aggregated across all samples).
    """
    function group_comparison_chart(dfs_per_group::Vector{DataFrame},
                                    scols_per_group::Vector{Vector{String}},
                                    group_names::Vector{String},
                                    output_path::String;
                                    top_n::Int = 15, rank = nothing,
                                    relative::Bool = true,
                                    subtitle = nothing,
                                    rank_order::Vector{String})
        isempty(dfs_per_group) && return nothing

        combined_df  = reduce((a, b) -> vcat(a, b; cols = :union), dfs_per_group)
        labels       = _label_rows(combined_df; rank, rank_order)
        unique_labels = String.(unique(labels))
        label_idx    = Dict(l => i for (i, l) in enumerate(unique_labels))
        n_labels     = length(unique_labels)
        n_groups     = length(group_names)

        counts = zeros(Float64, n_labels, n_groups)
        row_offset = 0
        for (gi, (df, scols)) in enumerate(zip(dfs_per_group, scols_per_group))
            for ri in 1:nrow(df)
                li = label_idx[labels[row_offset + ri]]
                for col in scols
                    v = df[ri, Symbol(col)]
                    counts[li, gi] += ismissing(v) ? 0.0 : Float64(v)
                end
            end
            row_offset += nrow(df)
        end

        totals = vec(sum(counts, dims = 2))
        order  = sortperm(totals, rev = true)
        if n_labels > top_n
            keep  = order[1:top_n]
            other = vec(sum(counts[order[(top_n + 1):end], :], dims = 1))
            counts       = vcat(counts[keep, :], other')
            final_labels = vcat(unique_labels[keep], ["Other"])
        else
            counts       = counts[order, :]
            final_labels = unique_labels[order]
        end

        if relative
            col_sums = vec(sum(counts, dims = 1))
            for j in 1:n_groups
                col_sums[j] > 0 && (counts[:, j] ./= col_sums[j])
            end
        end

        n_taxa = length(final_labels)
        colors = _palette_hex(n_taxa)
        _apply_grey!(colors, final_labels)

        ylabel = relative ? "Relative abundance" : "Read count"
        title  = relative ? "Group comparison" : "Group comparison (absolute)"

        traces = [Dict{String,Any}(
            "type"   => "bar",
            "name"   => final_labels[i],
            "x"      => group_names,
            "y"      => collect(counts[i, :]),
            "marker" => Dict("color" => colors[i]),
        ) for i in 1:n_taxa]

        layout = Dict{String,Any}(
            "barmode" => "stack",
            "title"   => Dict("text" => _title(title, subtitle)),
            "xaxis"   => Dict("title" => "Group", "tickangle" => -45),
            "yaxis"   => Dict{String,Any}("title" => ylabel),
        )
        relative && (layout["yaxis"]["range"] = [0, 1])

        mkpath(dirname(output_path))
        _write_json(output_path, Dict("data" => traces, "layout" => layout))
        return nothing
    end

end
