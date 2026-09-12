# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Route: /api/v1/studies/{study}/runs/{run}/analysis/read-funnel
#
# Per-sample read counts at every stage from raw reads to the final tables,
# with checks that steps which should keep every read did so.
using CSV, DataFrames

# Per-sample (processed, written) counts from cutadapt's stats file, keyed by the
# sample prefix of each block's -o output.
function _cutadapt_counts(stats_path::String)
    out = Dict{String,Tuple{Int,Int}}()
    isfile(stats_path) || return out
    num(s) = parse(Int, replace(s, "," => ""))
    for block in split(read(stats_path, String), "This is cutadapt")
        m_out = match(r"-o \S*?/([^/\s]+?)_R[12]_trimmed\.fastq\.gz", block)
        m_in  = match(r"Total (?:read pairs|reads) processed:\s+([\d,]+)", block)
        m_kept = match(r"(?:Pairs|Reads) written \(passing filters\):\s+([\d,]+)", block)
        (isnothing(m_out) || isnothing(m_in) || isnothing(m_kept)) && continue
        out[m_out.captures[1]] = (num(m_in.captures[1]), num(m_kept.captures[1]))
    end
    out
end

# Per-sample totals of a merged table, keyed by sample column.
function _table_sample_totals(csv_path::String, samples)
    isfile(csv_path) || return nothing
    df = CSV.read(csv_path, DataFrame)
    Dict(s => (s in names(df) ? sum(skipmissing(df[!, s]); init=0) : missing) for s in samples)
end

const _FUNNEL_STAGES = [
    ("raw",        "Raw reads"),
    ("trimmed",    "Primer-trimmed"),
    ("input",      "DADA2 input"),
    ("filtered",   "Quality-filtered"),
    ("denoisedF",  "Denoised (F)"),
    ("denoisedR",  "Denoised (R)"),
    ("merged",     "Pairs merged"),
    ("nochim",     "Chimera-free"),
    ("table",      "Merged table"),
    ("table_otu",  "OTU table"),
    ("table_cdhit", "CD-HIT table"),
]

function _read_funnel(run_dir::String)
    stats_csv = joinpath(run_dir, "dada2", "Tables", "pipeline_stats.csv")
    isfile(stats_csv) || return nothing
    stats = CSV.read(stats_csv, DataFrame)
    rename!(stats, names(stats)[1] => "sample")
    samples = String.(stats.sample)

    values = Dict{String,Dict{String,Any}}(s => Dict{String,Any}() for s in samples)
    for s in samples, c in names(stats)
        c == "sample" && continue
        v = stats[findfirst(==(s), samples), c]
        values[s][c] = ismissing(v) ? missing : Int(round(v))
    end
    cut = _cutadapt_counts(joinpath(run_dir, "cutadapt", "logs", "cutadapt_primer_trimming_stats.txt"))
    for s in samples
        haskey(cut, s) || continue
        values[s]["raw"], values[s]["trimmed"] = cut[s]
    end
    for (key, file) in (("table", "merged.csv"), ("table_otu", "merged_otu.csv"),
                        ("table_cdhit", "merged_cdhit.csv"))
        totals = _table_sample_totals(joinpath(run_dir, "merged", file), samples)
        isnothing(totals) && continue
        for s in samples
            values[s][key] = totals[s]
        end
    end

    present = [k for (k, _) in _FUNNEL_STAGES if any(s -> haskey(values[s], k), samples)]
    stages = [(; key=k, label=l) for (k, l) in _FUNNEL_STAGES if k in present]

    # Steps that move reads between files without filtering them must keep every read.
    checks = NamedTuple[]
    function equal_check(a, b, name)
        (a in present && b in present) || return
        bad = [s for s in samples
               if !ismissing(get(values[s], a, missing)) && !ismissing(get(values[s], b, missing)) &&
                  values[s][a] != values[s][b]]
        push!(checks, (; name, ok=isempty(bad),
                         detail=isempty(bad) ? "All $(length(samples)) samples match" :
                                "$(length(bad)) samples differ: $(join(first(bad, 5), ", "))$(length(bad) > 5 ? ", …" : "")"))
    end
    equal_check("trimmed", "input", "DADA2 reads every primer-trimmed read")
    equal_check("table", "table_otu", "OTU table holds the same reads as the merged table")
    equal_check("table", "table_cdhit", "CD-HIT table holds the same reads as the merged table")
    if "nochim" in present && "table" in present
        lost = sum(s -> coalesce(get(values[s], "nochim", 0), 0) - coalesce(get(values[s], "table", 0), 0), samples)
        push!(checks, (; name="Reads removed by merge_taxa filters", ok=true,
                         detail=lost == 0 ? "None" : "$lost reads across all samples"))
    end

    (; stages, samples=[(; sample=s, values=values[s]) for s in samples], checks)
end

@get "/api/v1/studies/{study}/runs/{run}/analysis/read-funnel" function(req, study::String, run::String)
    err = _validate_run_request(study, run)
    isnothing(err) || return err
    funnel = _read_funnel(_run_project_dir(study, run; group=_req_group(req)))
    isnothing(funnel) && return json_error(404, "no_stats", "No pipeline stats - run DADA2 first")
    json(funnel)
end
