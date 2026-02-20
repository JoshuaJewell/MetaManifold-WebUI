module TaxonomyTableTools

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

using CSV, DataFrames, Logging, YAML

export merge_taxonomy_counts, filter_table

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

    function build_taxonomy_table_rows(df::DataFrame)
        header = ["SeqName","Pident","Accession","rRNA","Organellum","specimen",
                "Domain","Supergroup","Division","Subdivision","Class","Order",
                "Family","Genus","Species"]
        out_rows = Vector{Vector{String}}(undef, nrow(df))

        for (i, r) in enumerate(eachrow(df))
            data = fill("", length(header))
            data[1] = string(r.SeqName)
            data[2] = string(r.Pident)

            tax_str = r.sseqtax
            if occursin('|', tax_str)
                # Full PR2 format: Accession|rRNA|Organellum|specimen|Domain|...|Species
                parts = split(tax_str, '|')
                for j in 1:min(length(parts), length(header) - 2)
                    data[j + 2] = parts[j]
                end
            else
                # Taxonomy-only format: Domain;Supergroup;...;Species[;]
                # No accession/rRNA/Organellum/specimen - those stay as empty strings
                parts = filter(!isempty, split(tax_str, ';'))
                for j in 1:min(length(parts), length(header) - 6)
                    data[j + 6] = parts[j]
                end
            end

            out_rows[i] = data
        end
        return header, out_rows
    end

    # Sorting helper
    function seqnum(x)
        if ismissing(x)
            return typemax(Int)
        end
        m = match(r"seq(\d+)", String(x))
        return m === nothing ? typemax(Int) : parse(Int, m.captures[1])
    end

    """
        merge_taxonomy_counts(taxonomy_vsearch_path, counts_csv_path)

    Reads taxonomy TSV from vsearch and counts CSV from DADA2 and returns a merged DataFrame.

    ## Arguments:
    - `taxonomy_vsearch_path` (): Path to vsearch output file.
    - `counts_csv_path`: Path to DADA2 idtax output file.
    """
    function merge_taxonomy_counts(
        taxonomy_vsearch_path,
        counts_csv_path)
        
        @info("Merging taxonomy counts.")

        # build taxonomy dataframe - FIXED VERSION
        imported = import_vsearch(taxonomy_vsearch_path)
        df_tax_raw = vsearch_to_df(imported)
        header, rows = build_taxonomy_table_rows(df_tax_raw)
        
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

        # Left join instead of outerjoin to keep all taxonomy rows
        merged_df = leftjoin(df_taxonomy, df_counts_prepared, on="SeqName")
        sort!(merged_df, "SeqName", by=seqnum)

        return merged_df
    end

    """
        filter_table(merged_df; filters_yaml_path, remove_empty_domain_override)

    Applies protist filtering to a merged DataFrame produced by merge_taxonomy_counts using rules defined in a YAML file.

    ## Arguments:
    - `merged_df`: Input DataFrame with taxonomy columns.
    - `filters_yaml_path` (default: "./inputs/protist_filter.yml"): Path to YAML file defining filters and mappings.
    - `remove_empty_domain_override` (default: nothing): Overrides YAML setting if provided.

    ## Expected YAML structure:
    ``` YAML
        mappings: { DivisionValue: SupergroupValue, ... }
        filters:
            - column: Domain
            pattern: Bacteria
        remove_empty_domain: true
    ```
    """
    function filter_table(
        merged_df::DataFrame,
        filters_yaml_path::String;
        remove_empty_domain_override::Union{Bool,Nothing} = nothing
    )
        # Load YAML config
        config = YAML.load_file(filters_yaml_path)

        df = deepcopy(merged_df)

        @info "Filtering table using configuration from $filters_yaml_path"

        # Load mappings
        mapping_config = get(config, "mappings", Dict())
        if !isempty(mapping_config)
            @assert mapping_config isa AbstractDict "mappings must be a dictionary in YAML"
        end

        mapping = Dict(string(k) => string(v) for (k,v) in mapping_config)

        # Apply Supergroup mapping if both columns exist
        if "Division" in names(df) && "Supergroup" in names(df) && !isempty(mapping)
            df.Supergroup = [
                ismissing(d) ? s : (haskey(mapping, string(d)) ? mapping[string(d)] : s)
                for (d, s) in zip(df.Division, df.Supergroup)
            ]
        elseif !isempty(mapping)
            missing_cols = filter(c -> c ∉ names(df), ["Division", "Supergroup"])
            @warn "Supergroup remapping skipped: column(s) not present in data: " *
                  "$(join(missing_cols, ", ")). The mappings block in your filter config " *
                  "may be tuned for a different reference database (e.g. PR2). " *
                  "Available columns: $(join(names(df), ", "))"
        end

        # Load remove_empty_domain flag
        remove_empty = get(config, "remove_empty_domain", true)
        if remove_empty_domain_override !== nothing
            remove_empty = remove_empty_domain_override
        end

        # Optionally remove blank/unassigned Domain entries
        if remove_empty && "Domain" in names(df)
            df = filter(row ->
                !ismissing(row.Domain) &&
                !isempty(strip(string(row.Domain))) &&
                lowercase(strip(string(row.Domain))) != "blank",
                df)
        elseif remove_empty
            @warn "remove_empty_domain skipped: 'Domain' column not present in data. " *
                  "This setting may be tuned for a different reference database (e.g. PR2). " *
                  "Available columns: $(join(names(df), ", "))"
        end

        # Load and apply filters
        raw_filters = get(config, "filters", [])
        @assert raw_filters isa Vector "filters must be a list in YAML"

        for item in raw_filters
            @assert item isa Dict "Each filter must be a key-value map"
            colname = get(item, "column", nothing)
            pattern = get(item, "pattern", nothing)
            @assert typeof(colname) <: AbstractString "filter column must be a string"
            @assert typeof(pattern) <: AbstractString "filter pattern must be a string"

            if !(colname in names(df))
                @warn "Skipping filter: column '$colname' not present in data. " *
                      "This filter may be tuned for a different reference database (e.g. PR2). " *
                      "Available columns: $(join(names(df), ", "))"
                continue
            end

            # Exclude rows where column contains the pattern (case-sensitive substring match)
            df = filter(row -> ismissing(row[colname]) || !occursin(pattern, string(row[colname])), df)
        end

        @info "Filtering complete. $(nrow(df))/$(nrow(merged_df)) rows retained."

        return df
    end

end