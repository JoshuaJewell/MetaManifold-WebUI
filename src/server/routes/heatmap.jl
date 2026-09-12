# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Heatmap palette shared by the results-table export and the publication tables.

using XLSX, DataFrames

const HEATMAP_MODES = ("none", "column", "table")
const _HEAT_WHITE = (255, 255, 255)
const _HEAT_HIGH  = (0x6F, 0xA8, 0xDC)   # light enough for black text at full strength
const _HEAT_LOW   = (0xE6, 0x84, 0x6B)   # negative end of a scale centred on zero

_heat_hex(rgb) = "#" * join(string(round(Int, c); base=16, pad=2) for c in rgb)
_heat_blend(to, t) = _heat_hex(Tuple(w + (c - w) * t for (w, c) in zip(_HEAT_WHITE, to)))

"""
    _heat_colour(v, m; diverging=false) -> Union{String,Nothing}

Hex fill for `v` on a scale from 0 to `m`. Zero, non-numeric values and empty
scales get no fill. With `diverging`, negative values use the low colour.
"""
function _heat_colour(v, m; diverging::Bool=false)
    (v isa Real && !(v isa Bool) && v != 0 && m > 0) || return nothing
    _heat_blend(diverging && v < 0 ? _HEAT_LOW : _HEAT_HIGH, min(abs(v) / m, 1.0))
end

_xlsx_fill!(sh, r::Integer, c::Integer, hex::AbstractString) =
    XLSX.setFill(sh, r, c; pattern="solid", fgColor="FF" * uppercase(hex[2:end]))

"""
    _count_fills(df, count_cols, mode) -> Vector{Tuple{Int,Int,String}}

(row, column, hex) fills for the count cells of `df`, row and column being
indices into `df`. `mode` "column" grades each count column from 0 to its own
maximum; "table" uses one maximum across all count columns.
"""
function _count_fills(df::AbstractDataFrame, count_cols::Vector{String}, mode::AbstractString)
    fills = Tuple{Int,Int,String}[]
    (mode == "none" || isempty(count_cols) || nrow(df) == 0) && return fills
    colmax(c) = maximum((abs(v) for v in df[!, c] if v isa Real); init=0.0)
    maxima = Dict(c => colmax(c) for c in count_cols)
    whole = maximum(values(maxima); init=0.0)
    for c in count_cols
        j = columnindex(df, c)
        m = mode == "table" ? whole : maxima[c]
        for (i, v) in enumerate(df[!, c])
            hex = _heat_colour(v, m)
            isnothing(hex) || push!(fills, (i, j, hex))
        end
    end
    fills
end
