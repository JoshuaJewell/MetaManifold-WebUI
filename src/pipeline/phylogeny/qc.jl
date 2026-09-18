# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## Quality summaries of each step
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

    # Support values are the numbers just after a closing bracket; IQ-TREE writes
    # "a/b" for two kinds of support, of which the first is kept.
    function _supports(newick::AbstractString)
        [parse(Float64, m.captures[1]) for m in eachmatch(r"\)([0-9]+(?:\.[0-9]+)?)(?:/[0-9.]+)?(?=[:,;)])", newick)]
    end

    function tree_qc(files)
        report = isfile(files["reference.iqtree"]) ? read(files["reference.iqtree"], String) : ""
        grab(re) = (m = match(re, report); isnothing(m) ? nothing : m.captures[1])
        num(re) = (v = grab(re); isnothing(v) ? nothing : parse(Float64, v))
        best = grab(r"Best-fit model according to \w+:\s*(\S+)")
        Dict{String,Any}(
            "kind"              => "tree",
            "model"             => something(best, grab(r"Model of substitution:\s*(\S+)"), Some(nothing)),
            "model_selected"    => !isnothing(best),
            "log_likelihood"    => num(r"Log-likelihood of the tree:\s*(-?[0-9.]+)"),
            "sequences"         => num(r"Input data:\s*(\d+) sequences"),
            "sites"             => num(r"Input data:\s*\d+ sequences with (\d+)"),
            "informative_sites" => num(r"Number of parsimony informative sites:\s*(\d+)"),
            "constant_sites"    => num(r"Number of constant sites:\s*(\d+)"),
            "supports"          => _supports(read(files["reference.treefile"], String)),
        )
    end

    function _jplace_queries(path)
        doc = JSON3.read(read(path, String))
        fields = String.(collect(doc.fields))
        lwr  = findfirst(==("like_weight_ratio"), fields)
        edge = findfirst(==("edge_num"), fields)
        out = Dict{String,Any}()
        for pl in doc.placements
            names = haskey(pl, :n) ? String.(collect(pl.n)) : [String(x[1]) for x in pl.nm]
            rows = collect(pl.p)
            best = isempty(rows) ? nothing : rows[argmax([Float64(r[lwr]) for r in rows])]
            for n in names
                out[n] = Dict{String,Any}(
                    "placements" => length(rows),
                    "best_lwr"   => isnothing(best) ? nothing : _round(Float64(best[lwr])),
                    "edge"       => isnothing(best) ? nothing : Int(best[edge]))
            end
        end
        out
    end

    function placement_qc(files)
        placed = _jplace_queries(files["placement.jplace"])
        raw = Dict(read_fasta(files["queries.fasta"]))
        trimmed = Dict(read_aligned(files["combined.trim.fasta"]))
        queries = map(sort(collect(keys(raw)))) do q
            t = get(trimmed, q, nothing)
            p = get(placed, q, nothing)
            Dict{String,Any}("name" => q, "residues" => length(raw[q]),
                             "trimmed_residues" => isnothing(t) ? 0 : count(!_is_gap, t),
                             "placed" => !isnothing(p),
                             "placements" => isnothing(p) ? 0 : p["placements"],
                             "best_lwr" => isnothing(p) ? nothing : p["best_lwr"])
        end
        Dict{String,Any}("kind" => "placement", "queries" => queries,
                         "placed" => count(q -> q["placed"], queries))
    end

    function accumulate_qc(files)
        before = _jplace_queries(files["placement.jplace"])
        after  = _jplace_queries(files["accumulated.jplace"])
        Dict{String,Any}("kind" => "accumulate",
                         "kept" => length(after),
                         "dropped" => sort([q for q in keys(before) if !haskey(after, q)]))
    end

    function step_qc(steps::Vector{Step}, dir, files, step::Step)
        placement = steps === PLACEMENT_STEPS
        queries = placement && isfile(files["queries.fasta"]) ?
                  Set(first.(read_fasta(files["queries.fasta"]))) : Set{String}()
        if step.name == "align"
            alignment_qc(files[step.outputs[1]]; queries)
        elseif step.name == "trim"
            src = placement ? "combined.aln.fasta" : "reference.aln.fasta"
            alignment_qc(files[src]; trimmed=files[step.outputs[1]],
                         columns=read_columns(files[step.outputs[2]]), queries)
        elseif step.name == "tree"
            tree_qc(files)
        elseif step.name == "place"
            placement_qc(files)
        else
            accumulate_qc(files)
        end
    end

    function _write_qc(steps, dir, files, step::Step)
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
