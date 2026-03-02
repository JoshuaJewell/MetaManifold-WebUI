# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

    function _generate_taxa_charts(df::DataFrame, scols::Vector{String},
                                    plots_dir::String, top_n::Int,
                                    subtitle::Union{Nothing,String},
                                    method_src_key::String;
                                    ranks, rank_order::Vector{String})
        mkpath(plots_dir)
        for (rankdir, rank) in ranks
            stem = "taxa_bar_$(method_src_key)_$(rankdir)"
            PipelinePlotsPlotly.taxa_bar_chart(df, scols,
                joinpath(plots_dir, "$(stem).json");
                top_n, rank, relative=true, subtitle, rank_order)
            PipelinePlotsPlotly.taxa_bar_chart(df, scols,
                joinpath(plots_dir, "$(stem)_absolute.json");
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
            print(io, rpad(PipelinePlotsPlotly._display_name(String(row.sample)), 20))
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
            print(io, rpad(PipelinePlotsPlotly._display_name(String(row.sample)), 20))
            print(io, rpad(string(row.richness), 12))
            print(io, rpad(string(round(row.shannon; digits=3)), 12))
            println(io, round(row.simpson; digits=3))
        end
    end
