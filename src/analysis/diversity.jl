module DiversityMetrics

# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

export richness, shannon, simpson

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

    Simpson diversity index 1 - sum(p_i^2), where p_i is the relative
    abundance of feature i.  Returns 0.0 when `counts` sums to zero.
    """
    function simpson(counts)
        n = sum(counts)
        n == 0 && return 0.0
        p = counts[counts .> 0] ./ n
        return 1.0 - sum(p .^ 2)
    end

end
