# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/studies/{study}/runs

using JSON3, Dates

const STAGES = [
    "cutadapt", "dada2_denoise", "dada2_classify", "swarm",
    "cdhit", "vsearch_asv", "vsearch_otu", "merge_taxa", "analyse_run"
]

# Maps each stage to the files whose existence (and mtime) indicate completion.
const STAGE_OUTPUTS = Dict(
    "cutadapt"        => run_dir -> joinpath(run_dir, "cutadapt"),
    "dada2_denoise"   => run_dir -> joinpath(run_dir, "dada2", "Checkpoints", "ckpt_chimera.RData"),
    "dada2_classify"  => run_dir -> joinpath(run_dir, "dada2", "Checkpoints", "checkpoint.RData"),
    "swarm"           => run_dir -> joinpath(run_dir, "swarm", "otu_table.csv"),
    "cdhit"           => run_dir -> joinpath(run_dir, "dada2", "Tables", "asvs.fasta"),
    "vsearch_asv"     => run_dir -> joinpath(run_dir, "vsearch", "asv_hits.tsv"),
    "vsearch_otu"     => run_dir -> joinpath(run_dir, "swarm", "vsearch", "otu_hits.tsv"),
    "merge_taxa"      => run_dir -> joinpath(run_dir, "merged", "merged.csv"),
    "analyse_run"     => run_dir -> joinpath(run_dir, "analysis", "pipeline_summary.csv"),
)

function _stage_status(stage::String, run_dir::String)
    sentinel = STAGE_OUTPUTS[stage](run_dir)
    if isdir(sentinel) || isfile(sentinel)
        # TODO: compare config.hash to detect stale outputs
        return "complete"
    end
    # Check if a job is currently running this stage
    running_jobs = filter(JobQueue.list_jobs()) do j
        j.stage == stage &&
        occursin(basename(run_dir), something(j.run, "")) &&
        j.status in (JobQueue.queued, JobQueue.running)
    end
    isempty(running_jobs) ? "not_started" : "running"
end

function _last_run(stage::String, run_dir::String)
    sentinel = STAGE_OUTPUTS[stage](run_dir)
    path = isdir(sentinel) ? sentinel : isfile(sentinel) ? sentinel : nothing
    isnothing(path) && return nothing
    string(unix2datetime(mtime(path))) * "Z"
end

function _run_json(study::String, run::String)
    run_dir     = joinpath(ServerState.projects_dir(), study, run)
    data_run    = joinpath(ServerState.data_dir(), study, run)
    samples     = _sample_names(study, run)

    stages = Dict(s => (;
        status   = _stage_status(s, run_dir),
        last_run = _last_run(s, run_dir),
    ) for s in STAGES)

    JSON3.write((;
        name    = run,
        study,
        sample_count = length(samples),
        samples,
        stages,
    ))
end

@get "/api/v1/studies/{study}/runs" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    runs = map(_run_names(study)) do run
        run_dir  = joinpath(ServerState.projects_dir(), study, run)
        samples  = _sample_names(study, run)
        stages   = Dict(s => (; status=_stage_status(s, run_dir)) for s in STAGES)
        (; name=run, study, sample_count=length(samples), stages)
    end
    json(runs)
end

@get "/api/v1/studies/{study}/runs/{run}" function(req, study::String, run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found in study '$study'")
    json(JSON3.read(_run_json(study, run)))
end
