# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

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
