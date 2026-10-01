# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## Quality summaries of the align and trim steps
# Written to qc/<step>.json after the step runs and read by the page.

    _is_gap(c::Char) = c == '-' || c == '.'
    _round(x) = round(x; digits=4)

    function read_columns(path::AbstractString)
        isfile(path) || return nothing
        m = match(r"#ColumnsMap\s+([0-9,\s]*)", read(path, String))
        isnothing(m) && return nothing
        [parse(Int, strip(x)) for x in split(m.captures[1], ',') if !isempty(strip(x))]
    end

    """
        alignment_qc(path; trimmed, columns, queries) -> Dict

    Column occupancy and per-sequence coverage of an alignment. With `columns`
    (trimAl's kept columns, from 0) and the `trimmed` alignment, it also says
    what trimming kept of each sequence and which sequences it removed.
    Occupancy is given for references and queries apart when `queries` names any.
    """
    function alignment_qc(path::AbstractString; trimmed=nothing, columns=nothing,
                          queries::AbstractSet=Set{String}())
        recs = read_aligned(path)
        L = isempty(recs) ? 0 : maximum(length(last(r)) for r in recs)
        kept = isnothing(columns) ? nothing : falses(L)
        isnothing(kept) || foreach(c -> c + 1 <= L && (kept[c + 1] = true), columns)
        occ(sel) = begin
            n = count(sel, recs)
            n == 0 && return Float64[]
            v = zeros(Int, L)
            for r in recs
                sel(r) || continue
                for (i, c) in enumerate(last(r))
                    _is_gap(c) || (v[i] += 1)
                end
            end
            _round.(v ./ n)
        end
        per = map(recs) do (name, seq)
            res = findall(!_is_gap, collect(seq))
            kept_res = isnothing(kept) ? nothing : count(i -> kept[i], res)
            Dict{String,Any}("name" => name, "query" => name in queries,
                             "residues" => length(res),
                             "kept_residues" => kept_res,
                             "span" => isempty(res) ? nothing : [first(res) - 1, last(res) - 1])
        end
        gaps = sum((count(_is_gap, last(r)) for r in recs); init=0)
        out = Dict{String,Any}(
            "kind"         => "alignment",
            "sequences"    => length(recs),
            "columns"      => L,
            "gap_fraction" => L == 0 ? 0.0 : _round(gaps / (L * length(recs))),
            "occupancy"    => occ(r -> true),
            "per_sequence" => per,
        )
        if !isempty(queries)
            out["occupancy_references"] = occ(r -> !(first(r) in queries))
            out["occupancy_queries"]    = occ(r -> first(r) in queries)
        end
        if !isnothing(columns)
            out["kept_columns"] = columns
            left = isnothing(trimmed) || !isfile(trimmed) ? nothing : Set(first.(read_aligned(trimmed)))
            out["removed"] = isnothing(left) ? String[] : [n for (n, _) in recs if !(n in left)]
        end
        out
    end

    const QC_STEPS = ("align", "trim")

    function step_qc(steps::Vector{Step}, dir, files, step::Step)
        placement = steps === PLACEMENT_STEPS
        queries = placement && isfile(files["queries.fasta"]) ?
                  Set(first.(read_fasta(files["queries.fasta"]))) : Set{String}()
        step.name == "align" && return alignment_qc(files[step.outputs[1]]; queries)
        src = placement ? "combined.aln.fasta" : "reference.aln.fasta"
        alignment_qc(files[src]; trimmed=files[step.outputs[1]],
                     columns=read_columns(files[step.outputs[2]]), queries)
    end

    function _write_qc(steps, dir, files, step::Step)
        step.name in QC_STEPS || return
        try
            qc = step_qc(steps, dir, files, step)
            mkpath(joinpath(dir, "qc"))
            write(joinpath(dir, "qc", "$(step.name).json"), JSON3.write(qc))
        catch e
            @warn "Phylogeny: no QC for $(step.name)" exception=e
        end
    end

    """
        trim_preview(steps, dir, files, trim) -> Dict

    Trim the step's input alignment with `trim` settings in a scratch directory
    and return its QC; nothing in `dir` changes.
    """
    function trim_preview(steps::Vector{Step}, dir, files, trim::AbstractDict)
        placement = steps === PLACEMENT_STEPS
        src = files[placement ? "combined.aln.fasta" : "reference.aln.fasta"]
        isfile(src) || error("Run the align step first")
        errors = Validation.ValidationError[]
        Validation._validate_trim(errors, _as_dict(trim), "trim", "trim")
        isempty(errors) || error(join((e.message for e in errors), "; "))
        mktempdir() do tmp
            out, cols = joinpath(tmp, "trim.fasta"), joinpath(tmp, "columns.txt")
            cmd = _words(tool_bin("trimal"), "-in", _sq(src), "-out", _sq(out), "-fasta",
                         trim_args(_as_dict(trim)), "-colnumbering", ">", _sq(cols))
            err = IOBuffer()
            ok = success(pipeline(`bash -lc $cmd`; stderr=err))
            ok || error("trimAl failed: $(strip(String(take!(err))))")
            queries = placement ? Set(first.(read_fasta(files["queries.fasta"]))) : Set{String}()
            alignment_qc(src; trimmed=out, columns=read_columns(cols), queries)
        end
    end
