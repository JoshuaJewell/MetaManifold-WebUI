module DiversityMetrics

# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

import Random

export richness, shannon, simpson, rarefy, srs, normalise_counts,
       auto_min_depth, alpha_diversity, NORMALISATION_METHODS,
       hellinger, transform_counts, TRANSFORM_METHODS

## Supported count depth-normalisation methods
const NORMALISATION_METHODS = ("none", "rarefy", "srs")

## Supported pre-dissimilarity transforms
const TRANSFORM_METHODS = ("none", "hellinger")

    """
        richness(counts) -> Int

    Observed richness: the number of non-zero features in `counts`.
    """
    richness(counts) = count(!iszero, counts)

    """
        shannon(counts) -> Float64

    Shannon diversity index H = -sum(p_i * ln(p_i)), where p_i is the
    relative abundance of feature i.  Zero-count features are ignored.
    Returns 0.0 when `counts` sums to zero.
    """
    function shannon(counts)
        n = sum(counts)
        n == 0 && return 0.0
        p = counts[counts .> 0] ./ n
        return -sum(p .* log.(p))
    end

    """
        simpson(counts) -> Float64

    Gini-Simpson diversity index 1 - sum(p_i^2), where p_i is the relative
    abundance of feature i.  Returns 0.0 when `counts` sums to zero.
    """
    function simpson(counts)
        n = sum(counts)
        n == 0 && return 0.0
        p = counts[counts .> 0] ./ n
        return 1.0 - sum(p .^ 2)
    end

    """
        rarefy(mat; depth, seed) -> Matrix{Float64}

    Subsample each sample row to exactly `depth` reads without replacement.
    Fractional input counts are rounded to the nearest integer before subsampling.
    """
    function rarefy(mat::Matrix{<:Real}; depth::Int, seed::Int)::Matrix{Float64}
        rng = Random.MersenneTwister(seed)
        nrows, nfeat = size(mat)
        out = zeros(Float64, nrows, nfeat)
        # One buffer reused across rows; grows to the largest library only.
        pool = Int[]
        for i in 1:nrows
            total = 0
            for j in 1:nfeat
                total += round(Int, mat[i, j])
            end
            resize!(pool, total)
            idx = 1
            for j in 1:nfeat
                for _ in 1:round(Int, mat[i, j])
                    pool[idx] = j
                    idx += 1
                end
            end
            # Partial Fisher-Yates: draw exactly `depth` reads without
            # replacement, tallying each as it is selected. Equivalent in
            # distribution to a full shuffle then taking the first `depth`,
            # but O(depth) rather than O(library size) random work.
            for k in 1:depth
                s = rand(rng, k:total)
                pool[k], pool[s] = pool[s], pool[k]
                out[i, pool[k]] += 1.0
            end
        end
        out
    end

    """
        srs(mat; depth, seed) -> Matrix{Float64}

    Scaling with Ranked Subsampling (Beule & Karlovsky 2020): normalise each
    sample row to exactly `depth` reads by scaling rather than resampling.

    Each count is multiplied by `depth / library_size` and floored; the shortfall
    between the floored total and `depth` is then handed out one read at a time to
    the features with the largest discarded fractions.  Ties on the fraction are
    broken by the larger scaled integer part, and any tie surviving that is broken
    at random from `seed` — the only point at which SRS consults the RNG, which is
    why it reproduces far more of the original community structure than
    rarefaction at the same depth.

    Only features the sample actually observed can receive a remainder, so SRS
    never invents a read for an absent taxon.  Rows whose library size is zero are
    left as zeros.
    """
    function srs(mat::Matrix{<:Real}; depth::Int, seed::Int)::Matrix{Float64}
        rng = Random.MersenneTwister(seed)
        nrows, nfeat = size(mat)
        out = zeros(Float64, nrows, nfeat)
        # Scratch buffers reused across rows.
        parts = zeros(Int, nfeat)
        fracs = zeros(Float64, nfeat)
        jitter = zeros(Float64, nfeat)
        candidates = Int[]

        for i in 1:nrows
            total = 0.0
            for j in 1:nfeat
                total += mat[i, j]
            end
            total > 0 || continue

            scale = depth / total
            assigned = 0
            empty!(candidates)
            for j in 1:nfeat
                scaled = mat[i, j] * scale
                whole = floor(Int, scaled)
                parts[j] = whole
                fracs[j] = scaled - whole
                jitter[j] = rand(rng)
                out[i, j] = Float64(whole)
                assigned += whole
                mat[i, j] > 0 && push!(candidates, j)
            end

            deficit = depth - assigned
            deficit > 0 || continue
            sort!(candidates; by = j -> (-fracs[j], -parts[j], jitter[j]))
            for k in 1:min(deficit, length(candidates))
                out[i, candidates[k]] += 1.0
            end
        end
        out
    end

    """
        hellinger(mat) -> Matrix{Float64}

    Hellinger transform: the element-wise square root of each row's relative
    abundances.  It down-weights the dominant taxa that would otherwise drive
    the ordination.  Euclidean distance on the result is the Hellinger distance.
    Rows summing to zero are returned as zeros.
    """
    function hellinger(mat::Matrix{<:Real})::Matrix{Float64}
        out = zeros(Float64, size(mat))
        for i in axes(mat, 1)
            total = 0.0
            for j in axes(mat, 2)
                total += mat[i, j]
            end
            total > 0 || continue
            for j in axes(mat, 2)
                v = mat[i, j]
                out[i, j] = v > 0 ? sqrt(v / total) : 0.0
            end
        end
        out
    end

    """
        transform_counts(mat; method) -> Matrix{Float64}

    Dispatcher for the pre-dissimilarity transform. `method` is one of `"none"`
    or `"hellinger"`.
    """
    function transform_counts(mat::Matrix{<:Real}; method::String)::Matrix{Float64}
        method in TRANSFORM_METHODS || error(
            "Unknown transform method: $method " *
            "(expected one of $(join(TRANSFORM_METHODS, ", ")))")
        method == "none" ? Matrix{Float64}(mat) : hellinger(mat)
    end

    """
        auto_min_depth(lib_sizes) -> Int

    Auto rarefaction depth: the minimum strictly-positive library size across
    `lib_sizes`, or 0 when no sample has any reads.
    """
    function auto_min_depth(lib_sizes)
        best = nothing
        for s in lib_sizes
            s > 0 || continue
            (best === nothing || s < best) && (best = s)
        end
        best === nothing ? 0 : Int(best)
    end

    """
        normalise_counts(mat; method, depth, seed)
            -> (; mat::Matrix{Float64}, kept::Vector{Int})

    Dispatcher for count depth-normalisation. `method` is one of `"none"`,
    `"rarefy"` (random subsampling without replacement) or `"srs"` (scaling with
    ranked subsampling).

    `depth = 0` selects auto mode: the resolved depth is the minimum library
    size across samples that have at least one read.  Samples whose library
    size is strictly below the resolved depth are dropped before normalisation:
    neither method can scale a short library up to the target depth.

    Returns a named tuple `(; mat, kept)` where `kept` is the 1-based vector
    of retained row indices.  Callers must re-index any parallel label vectors
    (sample names, group labels, etc.) using `kept`.
    """
    function normalise_counts(mat::Matrix{<:Real};
                              method::String,
                              depth::Int,
                              seed::Int)::NamedTuple{(:mat, :kept), Tuple{Matrix{Float64}, Vector{Int}}}
        method in NORMALISATION_METHODS || error(
            "Unknown normalisation method: $method " *
            "(expected one of $(join(NORMALISATION_METHODS, ", ")))")
        lib_sizes = vec(sum(mat; dims=2))
        method == "none" && return (; mat=Matrix{Float64}(mat), kept=collect(1:size(mat, 1)))

        resolved_depth = depth == 0 ? auto_min_depth(lib_sizes) : depth

        @info "Normalisation: method=$method, resolved depth=$resolved_depth ($(length(lib_sizes)) samples, lib sizes $(Int.(extrema(lib_sizes))))"
        kept = findall(>=(resolved_depth), lib_sizes)
        dropped = size(mat, 1) - length(kept)
        dropped > 0 && @info "Normalisation: dropped $dropped samples below depth $resolved_depth"

        sub = mat[kept, :]
        normalised = method == "srs" ? srs(sub; depth=resolved_depth, seed) :
                                       rarefy(sub; depth=resolved_depth, seed)
        (; mat=normalised, kept)
    end

    """
        alpha_diversity(mat; method, depth, seed, iterations=1)
            -> (; kept, richness, shannon, simpson)

    Richness, Shannon and Simpson per retained sample. With rarefaction and
    `iterations > 1`, each metric is the mean over that many independent draws
    (seeds `seed`, `seed + 1`, ...), so one draw's chance does not decide it.
    """
    function alpha_diversity(mat::Matrix{<:Real}; method::String, depth::Int, seed::Int,
                             iterations::Int=1)
        draws = method == "rarefy" ? max(1, iterations) : 1
        kept = Int[]
        r = Float64[]; sh = Float64[]; si = Float64[]
        for k in 1:draws
            norm = normalise_counts(mat; method, depth, seed=seed + k - 1)
            if k == 1
                kept = norm.kept
                r  = zeros(length(kept)); sh = zeros(length(kept)); si = zeros(length(kept))
            end
            for i in eachindex(kept)
                counts = round.(Int, norm.mat[i, :])
                r[i]  += richness(counts) / draws
                sh[i] += shannon(counts) / draws
                si[i] += simpson(counts) / draws
            end
        end
        (; kept, richness=r, shannon=sh, simpson=si)
    end

end
