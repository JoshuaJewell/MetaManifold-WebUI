module PipelinePlots

# CairoMakie-based plotting for the analysis stage.
#
# All public functions write a PDF to the given path and return nothing.
# Module is named PipelinePlots (not Plots) to avoid shadowing Plots.jl.
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export taxa_bar_chart, filter_composition_plot, alpha_diversity_plot,
       nmds_plot, alpha_boxplot, pipeline_stats_plot, group_comparison_chart

    using CairoMakie, DataFrames

    CairoMakie.activate!(type = "pdf")

    # Colour palette
    # Return `n` visually distinct colours.  Starts from Makie's
    # Wong palette (7 colours optimised for colour-blindness) then
    # extends via golden-angle hue rotation.
    function _palette(n::Int)
        base = Makie.wong_colors()
        n <= length(base) && return collect(base[1:n])
        colors = collect(base)
        for i in (length(base) + 1):n
            hue = mod(i * 137.508, 360.0)
            # Simple HSV to RGB with S=0.65, V=0.85
            s, v = 0.65, 0.85
            c = v * s
            x = c * (1.0 - abs(mod(hue / 60.0, 2.0) - 1.0))
            m = v - c
            r, g, b = if hue < 60
                (c + m, x + m, m)
            elseif hue < 120
                (x + m, c + m, m)
            elseif hue < 180
                (m, c + m, x + m)
            elseif hue < 240
                (m, x + m, c + m)
            elseif hue < 300
                (x + m, m, c + m)
            else
                (c + m, m, x + m)
            end
            push!(colors, RGBAf(r, g, b, 1.0))
        end
        return colors[1:n]
    end

    ## Taxonomy label helpers
    # Return the value of the most specific (lowest) non-empty taxonomy rank.
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
        return label
    end

    # Return the value at a specified rank, or "Unclassified".
    function _rank_label(row, rank::String)
        sym = Symbol(rank)
        if hasproperty(row, sym)
            val = row[sym]
            if !ismissing(val) && !isempty(strip(string(val)))
                return string(val)
            end
        end
        return "Unclassified"
    end

    # Assign a taxon label to every row of `df`.
    function _label_rows(df::DataFrame; rank=nothing,
                         rank_order::Vector{String})::Vector{String}
        if isnothing(rank)
            String[_lowest_rank_label(row; rank_order) for row in eachrow(df)]
        else
            String[_rank_label(row, rank) for row in eachrow(df)]
        end
    end

    ## Display-name helper
    # Strip FASTQ suffixes so tick labels show sample stems only.
    const _FASTQ_SUFFIX = r"_R[12](?:_filt|_trimmed)?\.fastq\.gz$"
    _display_name(name::String) = replace(name, _FASTQ_SUFFIX => "")
    _display_names(names::Vector{String}) = String[_display_name(n) for n in names]

    ## Subtitle helper
    function _add_subtitle!(fig, subtitle)
        isnothing(subtitle) && return
        Label(fig[0, 1:end], subtitle; fontsize=11, color=:gray40)
    end

    # Overrides Unclassified / Other as grey
    const _GREY = RGBAf(0.7, 0.7, 0.7, 1.0)

    function _apply_grey_override!(colors::Vector, labels::Vector{String})
        for (i, l) in enumerate(labels)
            if l in ("Unclassified", "Other")
                colors[i] = _GREY
            end
        end
    end

    ## Stacked-bar data builder
    # From a label-per-row vector and sample columns, produce the
    # aggregated counts matrix (labels x samples), final label list,
    # and per-column normalisation (relative or absolute).
    function _aggregate_taxa(df::DataFrame, sample_cols::Vector{String},
                             labels::AbstractVector{<:AbstractString};
                             top_n::Int=15, relative::Bool=true)
        unique_labels = String.(unique(labels))
        label_idx = Dict(l => i for (i, l) in enumerate(unique_labels))
        n_labels  = length(unique_labels)
        n_samples = length(sample_cols)

        counts = zeros(Float64, n_labels, n_samples)
        for (ri, row) in enumerate(eachrow(df))
            li = label_idx[labels[ri]]
            for (j, col) in enumerate(sample_cols)
                v = row[Symbol(col)]
                counts[li, j] += ismissing(v) ? 0.0 : Float64(v isa AbstractString ? parse(Float64, v) : v)
            end
        end

        # Sort labels by total abundance (descending).
        totals = vec(sum(counts, dims=2))
        order  = sortperm(totals, rev=true)

        # Collapse beyond top_n into "Other".
        if n_labels > top_n
            keep  = order[1:top_n]
            other = vec(sum(counts[order[(top_n + 1):end], :], dims=1))
            counts       = vcat(counts[keep, :], other')
            final_labels = vcat(unique_labels[keep], ["Other"])
        else
            counts       = counts[order, :]
            final_labels = unique_labels[order]
        end

        # Optionally normalise to relative abundance.
        if relative
            col_sums = vec(sum(counts, dims=1))
            for j in 1:n_samples
                col_sums[j] > 0 && (counts[:, j] ./= col_sums[j])
            end
        end

        return counts, final_labels
    end

    # Build the three parallel vectors (x, y, grp) for a stacked barplot.
    function _long_format(counts::Matrix{Float64})
        n_grp, n_x = size(counts)
        x   = Int[]
        y   = Float64[]
        grp = Int[]
        for j in 1:n_x, i in 1:n_grp
            push!(x, j)
            push!(y, counts[i, j])
            push!(grp, i)
        end
        return x, y, grp
    end

    ## Public: taxa_bar_chart
    """
        taxa_bar_chart(df, sample_cols, output_pdf; top_n=15, rank=nothing,
                       relative=true, subtitle=nothing)

    Stacked bar chart of taxonomic composition per sample.

    `df` is a merged/filtered CSV read as a DataFrame.  `sample_cols` lists
    the column names that hold per-sample ASV counts.  Taxa are labelled by
    the lowest assigned taxonomic rank (or `rank` if specified) and the
    `top_n` most abundant are shown; the rest are collapsed to "Other".

    When `relative=false`, shows absolute read counts instead of proportions.
    """
    function taxa_bar_chart(df::DataFrame, sample_cols::Vector{String},
                            output_pdf::String; top_n::Int=15, rank=nothing,
                            relative::Bool=true,
                            subtitle::Union{Nothing,String}=nothing,
                            rank_order::Vector{String})
        (nrow(df) == 0 || isempty(sample_cols)) && return nothing

        labels = _label_rows(df; rank, rank_order)
        counts, final_labels = _aggregate_taxa(df, sample_cols, labels;
                                               top_n, relative)
        n_taxa    = length(final_labels)
        n_samples = length(sample_cols)
        x, y, grp = _long_format(counts)

        colors    = _palette(n_taxa)
        _apply_grey_override!(colors, final_labels)
        color_vec = [colors[g] for g in grp]

        display_names = _display_names(sample_cols)
        ylabel = relative ? "Relative abundance" : "Read count"
        title  = relative ? "Taxonomic composition" : "Taxonomic composition (absolute)"

        fig = Figure(size = (max(800, n_samples * 55 + 250), 600))
        ax  = Axis(fig[1, 1];
                    xlabel = "Sample",
                    ylabel = ylabel,
                    xticks = (1:n_samples, display_names),
                    xticklabelrotation = π / 4,
                    title  = title)

        barplot!(ax, x, y; stack=grp, color=color_vec)
        relative && ylims!(ax, 0, 1)

        elements = [PolyElement(color=colors[i]) for i in 1:n_taxa]
        Legend(fig[1, 2], elements, final_labels; framevisible=false,
               nbanks=1, labelsize=11)

        _add_subtitle!(fig, subtitle)
        save(output_pdf, fig)
        return nothing
    end

    ## Public: filter_composition_plot
    """
        filter_composition_plot(filter_totals, sample_names, output_pdf;
                                subtitle=nothing)

    Stacked bar chart showing absolute read counts per sample, one colour
    per biological-group filter (e.g. "protists", "bacteria", "Unclassified").

    `filter_totals` maps filter name to a Vector{Float64} of per-sample totals.
    """
    function filter_composition_plot(filter_totals::Dict{String, Vector{Float64}},
                                     sample_names::Vector{String},
                                     output_pdf::String;
                                     subtitle::Union{Nothing,String}=nothing,
                                     colour_overrides::Dict{String,String}=Dict{String,String}())
        isempty(filter_totals) && return nothing
        n_samples = length(sample_names)

        display_names = _display_names(sample_names)
        filter_names = sort(collect(keys(filter_totals)))
        # Move "Unclassified" to the end if present.
        if "Unclassified" in filter_names
            filter!(!=( "Unclassified"), filter_names)
            push!(filter_names, "Unclassified")
        end
        n_filters = length(filter_names)

        x   = Int[]
        y   = Float64[]
        grp = Int[]
        for j in 1:n_samples, (i, fname) in enumerate(filter_names)
            push!(x, j)
            push!(y, filter_totals[fname][j])
            push!(grp, i)
        end

        colors    = _palette(n_filters)
        _apply_grey_override!(colors, filter_names)
        # Apply per-filter colour overrides from YAML configs.
        for (i, fname) in enumerate(filter_names)
            if haskey(colour_overrides, fname)
                hex = colour_overrides[fname]
                try
                    colors[i] = Makie.RGBAf(Makie.Colors.parse(Makie.Colors.Colorant, hex))
                catch
                    @warn "Invalid colour '$hex' for filter '$fname', using default"
                end
            end
        end
        color_vec = [colors[g] for g in grp]

        fig = Figure(size = (max(800, n_samples * 55 + 250), 600))
        ax  = Axis(fig[1, 1];
                    xlabel = "Sample",
                    ylabel = "Read count",
                    xticks = (1:n_samples, display_names),
                    xticklabelrotation = π / 4,
                    title  = "Filter composition")

        barplot!(ax, x, y; stack=grp, color=color_vec)

        elements = [PolyElement(color=colors[i]) for i in 1:n_filters]
        Legend(fig[1, 2], elements, filter_names; framevisible=false)

        _add_subtitle!(fig, subtitle)
        save(output_pdf, fig)
        return nothing
    end

    ## Public: alpha_diversity_plot
    """
        alpha_diversity_plot(alpha_df, output_pdf; subtitle=nothing)

    Three-panel bar chart (richness, Shannon, Simpson) with samples on X.

    `alpha_df` must have columns `:sample`, `:richness`, `:shannon`, `:simpson`.
    """
    function alpha_diversity_plot(alpha_df::DataFrame, output_pdf::String;
                                  subtitle::Union{Nothing,String}=nothing)
        nrow(alpha_df) == 0 && return nothing
        n = nrow(alpha_df)
        samples = _display_names(String.(alpha_df.sample))

        metrics = [
            ("Richness (observed ASVs)", Float64.(alpha_df.richness)),
            ("Shannon index",            Float64.(alpha_df.shannon)),
            ("Simpson index",            Float64.(alpha_df.simpson)),
        ]

        fig = Figure(size = (max(700, n * 50), 900))
        for (row, (title, vals)) in enumerate(metrics)
            ax = Axis(fig[row, 1];
                      ylabel = title,
                      xticks = (1:n, samples),
                      xticklabelrotation = π / 4)
            row == 1 && (ax.title = "Alpha diversity")
            row < 3  && (ax.xticklabelsvisible = false; ax.xlabelvisible = false)
            row == 3 && (ax.xlabel = "Sample")
            barplot!(ax, 1:n, vals; color=Makie.wong_colors()[1])
        end

        _add_subtitle!(fig, subtitle)
        save(output_pdf, fig)
        return nothing
    end

    ## Public: nmds_plot
    # Marker shapes for dual-encoded NMDS.
    const _MARKER_SHAPES = [:circle, :rect, :diamond, :utriangle,
                            :dtriangle, :cross, :star5, :hexagon]

    """
        nmds_plot(coords, labels, output_pdf; colour_by=nothing,
                  shape_by=nothing, colour_label="Group",
                  shape_label="Group", stress=nothing, subtitle=nothing)

    NMDS ordination scatter plot.

    `coords` is an nx2 matrix (from R `metaMDS\$points`).
    `labels` are sample names.  When `colour_by` is a `Vector{String}`,
    points are coloured by group membership.  When `shape_by` is also
    provided, points are additionally encoded by marker shape with a
    grouped legend.
    """
    function nmds_plot(coords::Matrix{Float64}, labels::Vector{String},
                       output_pdf::String;
                       colour_by::Union{Nothing, Vector{String}}=nothing,
                       shape_by::Union{Nothing, Vector{String}}=nothing,
                       colour_label::String="Group",
                       shape_label::String="Group",
                       stress::Union{Nothing, Float64}=nothing,
                       subtitle::Union{Nothing,String}=nothing)
        n = size(coords, 1)
        n == 0 && return nothing

        fig = Figure(size = (700, 650))
        ax  = Axis(fig[1, 1];
                    xlabel = "NMDS1", ylabel = "NMDS2",
                    title  = "NMDS ordination")

        if isnothing(colour_by) && isnothing(shape_by)
            # Plain scatter - no grouping.
            scatter!(ax, coords[:, 1], coords[:, 2];
                     markersize=10, color=Makie.wong_colors()[1])

        elseif !isnothing(colour_by) && isnothing(shape_by)
            # Colour-only grouping.
            cgroups = unique(colour_by)
            colors  = _palette(length(cgroups))
            for (gi, g) in enumerate(cgroups)
                mask = colour_by .== g
                scatter!(ax, coords[mask, 1], coords[mask, 2];
                         markersize=10, color=colors[gi], label=g)
            end
            Legend(fig[1, 2], ax; framevisible=false)

        else
            # Dual-encoded: colour x shape.
            cgroups = isnothing(colour_by) ? ["all"] : unique(colour_by)
            sgroups = unique(shape_by)
            colors  = _palette(length(cgroups))
            cb      = isnothing(colour_by) ? fill("all", n) : colour_by

            # Plot each combination.
            for (ci, cg) in enumerate(cgroups), (si, sg) in enumerate(sgroups)
                mask = (cb .== cg) .& (shape_by .== sg)
                any(mask) || continue
                mi = mod1(si, length(_MARKER_SHAPES))
                scatter!(ax, coords[mask, 1], coords[mask, 2];
                         markersize=10, color=colors[ci],
                         marker=_MARKER_SHAPES[mi])
            end

            # Build grouped legend: elements, labels, titles.
            colour_elements = [MarkerElement(color=colors[i], marker=:circle,
                                             markersize=10)
                               for i in 1:length(cgroups)]
            colour_labels   = [string(g) for g in cgroups]

            shape_elements  = [MarkerElement(color=:gray50,
                                             marker=_MARKER_SHAPES[mod1(i, length(_MARKER_SHAPES))],
                                             markersize=10)
                               for i in 1:length(sgroups)]
            shape_labels    = [string(g) for g in sgroups]

            Legend(fig[1, 2],
                   [colour_elements, shape_elements],
                   [colour_labels, shape_labels],
                   [colour_label, shape_label];
                   framevisible=false)
        end

        # Annotate stress value.
        if !isnothing(stress)
            text!(ax, 0.02, 0.98;
                  text  = "stress = $(round(stress; digits=3))",
                  space = :relative,
                  align = (:left, :top),
                  fontsize = 12)
        end

        _add_subtitle!(fig, subtitle)
        save(output_pdf, fig)
        return nothing
    end

    ## Public: alpha_boxplot (group / study level)
    """
        alpha_boxplot(all_alpha, group_labels, output_pdf; subtitle=nothing)

    Three-panel boxplot (richness, Shannon, Simpson) comparing groups.

    `all_alpha` is a `Vector{DataFrame}`, one per group, each with columns
    `:sample`, `:richness`, `:shannon`, `:simpson`.
    `group_labels` names each group.
    """
    function alpha_boxplot(all_alpha::Vector{DataFrame},
                           group_labels::Vector{String},
                           output_pdf::String;
                           subtitle::Union{Nothing,String}=nothing)
        isempty(all_alpha) && return nothing
        n_groups = length(group_labels)
        colors   = _palette(n_groups)

        metrics = [:richness, :shannon, :simpson]
        titles  = ["Richness (observed ASVs)", "Shannon index", "Simpson index"]

        fig = Figure(size = (max(600, n_groups * 120), 900))
        for (row, (col, title)) in enumerate(zip(metrics, titles))
            ax = Axis(fig[row, 1];
                      ylabel = title,
                      xticks = (1:n_groups, group_labels),
                      xticklabelrotation = π / 6)
            row == 1 && (ax.title = "Alpha diversity comparison")
            row < 3  && (ax.xticklabelsvisible = false; ax.xlabelvisible = false)
            row == 3 && (ax.xlabel = "Group")

            positions = Int[]
            values    = Float64[]
            cols_vec  = RGBAf[]
            for (gi, adf) in enumerate(all_alpha)
                vals = Float64.(adf[!, col])
                append!(positions, fill(gi, length(vals)))
                append!(values, vals)
                append!(cols_vec, fill(colors[gi], length(vals)))
            end

            boxplot!(ax, positions, values; color=cols_vec)
        end

        _add_subtitle!(fig, subtitle)
        save(output_pdf, fig)
        return nothing
    end

    ## Public: pipeline_stats_plot
    # Readable stage labels for pipeline_stats columns.
    const _STAGE_LABELS = Dict(
        "input"     => "Input",
        "filtered"  => "Filtered",
        "denoisedF" => "Denoised (F)",
        "denoisedR" => "Denoised (R)",
        "merged"    => "Merged",
        "nochim"    => "Chimera-free",
    )

    # Canonical stage order - extras (reads_*) are appended dynamically.
    const _STAGE_ORDER = ["input", "filtered", "denoisedF", "denoisedR",
                          "merged", "nochim"]

    function _stage_label(col::String)
        haskey(_STAGE_LABELS, col) && return _STAGE_LABELS[col]
        # reads_protist_filter -> "After protist_filter"
        if startswith(col, "reads_")
            return "After " * replace(col[7:end], "_" => " ")
        end
        return col
    end

    """
        pipeline_stats_plot(stats_df, output_pdf; subtitle=nothing)

    Grouped bar chart showing read counts at each pipeline stage per sample.

    `stats_df` must have a `:sample` column and one or more stage columns
    (input, filtered, denoisedF, ..., reads_*).
    """
    function pipeline_stats_plot(stats_df::DataFrame, output_pdf::String;
                                  subtitle::Union{Nothing,String}=nothing)
        nrow(stats_df) == 0 && return nothing

        all_cols = names(stats_df)

        # Select stage columns in canonical order, then reads_* sorted.
        stage_cols = String[]
        for col in _STAGE_ORDER
            col in all_cols && push!(stage_cols, col)
        end
        reads_cols = sort([c for c in all_cols
                           if startswith(c, "reads_")])
        append!(stage_cols, reads_cols)
        isempty(stage_cols) && return nothing

        sample_names = _display_names(String.(stats_df.sample))
        n_samples = length(sample_names)
        n_stages  = length(stage_cols)
        stage_labels = [_stage_label(c) for c in stage_cols]

        colors = _palette(n_samples)

        # Build dodge-grouped barplot data.
        xs   = Int[]
        ys   = Float64[]
        grps = Int[]
        for (si, row) in enumerate(eachrow(stats_df))
            for (xi, col) in enumerate(stage_cols)
                v = row[Symbol(col)]
                push!(xs, xi)
                push!(ys, ismissing(v) ? 0.0 : Float64(v))
                push!(grps, si)
            end
        end

        fig = Figure(size = (max(800, n_stages * 100 + 200), 600))
        ax  = Axis(fig[1, 1];
                    xlabel = "Pipeline stage",
                    ylabel = "Read count",
                    xticks = (1:n_stages, stage_labels),
                    xticklabelrotation = π / 4,
                    title  = "Reads through pipeline stages")

        barplot!(ax, xs, ys; dodge=grps,
                 color=[colors[g] for g in grps])

        elements = [PolyElement(color=colors[i]) for i in 1:n_samples]
        Legend(fig[1, 2], elements, sample_names; framevisible=false,
               labelsize=11)

        _add_subtitle!(fig, subtitle)
        save(output_pdf, fig)
        return nothing
    end

    ## Public: group_comparison_chart
    """
        group_comparison_chart(dfs_per_group, scols_per_group, group_names,
                               output_pdf; top_n=15, rank=nothing,
                               relative=true, subtitle=nothing)

    Stacked bar chart with one bar per group (aggregated reads across all
    samples in the group).  Useful for cross-group comparison at a glance.
    """
    function group_comparison_chart(dfs_per_group::Vector{DataFrame},
                                    scols_per_group::Vector{Vector{String}},
                                    group_names::Vector{String},
                                    output_pdf::String;
                                    top_n::Int=15, rank=nothing,
                                    relative::Bool=true,
                                    subtitle::Union{Nothing,String}=nothing,
                                    rank_order::Vector{String})
        isempty(dfs_per_group) && return nothing

        # Combine all DataFrames and build labels.
        combined_df = reduce((a, b) -> vcat(a, b; cols=:union), dfs_per_group)
        labels = _label_rows(combined_df; rank, rank_order)

        # Build a synthetic count matrix: one column per group.
        unique_labels = String.(unique(labels))
        label_idx = Dict(l => i for (i, l) in enumerate(unique_labels))
        n_labels = length(unique_labels)
        n_groups = length(group_names)

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

        # Sort by total abundance, collapse beyond top_n.
        totals = vec(sum(counts, dims=2))
        order = sortperm(totals, rev=true)
        if n_labels > top_n
            keep = order[1:top_n]
            other = vec(sum(counts[order[(top_n + 1):end], :], dims=1))
            counts = vcat(counts[keep, :], other')
            final_labels = vcat(unique_labels[keep], ["Other"])
        else
            counts = counts[order, :]
            final_labels = unique_labels[order]
        end

        if relative
            col_sums = vec(sum(counts, dims=1))
            for j in 1:n_groups
                col_sums[j] > 0 && (counts[:, j] ./= col_sums[j])
            end
        end

        n_taxa = length(final_labels)
        x, y, grp = _long_format(counts)

        colors = _palette(n_taxa)
        _apply_grey_override!(colors, final_labels)
        color_vec = [colors[g] for g in grp]

        ylabel = relative ? "Relative abundance" : "Read count"
        title = relative ? "Group comparison" : "Group comparison (absolute)"

        fig = Figure(size = (max(600, n_groups * 80 + 250), 600))
        ax = Axis(fig[1, 1];
                   xlabel = "Group",
                   ylabel = ylabel,
                   xticks = (1:n_groups, group_names),
                   xticklabelrotation = π / 4,
                   title = title)

        barplot!(ax, x, y; stack=grp, color=color_vec)
        relative && ylims!(ax, 0, 1)

        elements = [PolyElement(color=colors[i]) for i in 1:n_taxa]
        Legend(fig[1, 2], elements, final_labels; framevisible=false,
               nbanks=1, labelsize=11)

        _add_subtitle!(fig, subtitle)
        save(output_pdf, fig)
        return nothing
    end

end
