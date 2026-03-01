# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/studies  and  /api/v1/init

using JSON3

## helpers

function _study_names()
    isdir(ServerState.data_dir()) || return String[]
    filter(readdir(ServerState.data_dir())) do name
        isdir(joinpath(ServerState.data_dir(), name))
    end
end

function _run_names(study::String)
    d = joinpath(ServerState.data_dir(), study)
    isdir(d) || return String[]
    filter(readdir(d)) do name
        isdir(joinpath(d, name)) &&
        !isempty(filter(f -> endswith(f, ".fastq.gz"), readdir(joinpath(d, name))))
    end
end

function _group_names(study::String)
    d = joinpath(ServerState.data_dir(), study)
    isdir(d) || return String[]
    # A group directory contains subdirectories (runs), not fastq files directly.
    filter(readdir(d)) do name
        sub = joinpath(d, name)
        isdir(sub) &&
        isempty(filter(f -> endswith(f, ".fastq.gz"), readdir(sub))) &&
        !isempty(filter(n -> isdir(joinpath(sub, n)), readdir(sub)))
    end
end

function _sample_names(study::String, run::String)
    d = joinpath(ServerState.data_dir(), study, run)
    isdir(d) || return String[]
    files = filter(f -> endswith(f, "_R1.fastq.gz") || endswith(f, "_R1_001.fastq.gz"),
                   readdir(d))
    map(f -> replace(f, r"_R1(_001)?\.fastq\.gz$" => ""), files)
end

function _active_job_count(study::String)
    jobs = JobQueue.list_jobs(; study)
    count(j -> j.status in (JobQueue.queued, JobQueue.running), jobs)
end

function _study_json(name::String)
    runs   = _run_names(name)
    groups = _group_names(name)
    JSON3.write((;
        name,
        run_count        = length(runs),
        group_count      = length(groups),
        active_job_count = _active_job_count(name),
        runs,
        groups,
    ))
end

function _ensure_all_projects()
    initialised = String[]
    for study in _study_names()
        proj_dir = joinpath(ServerState.projects_dir(), study)
        isdir(proj_dir) && continue
        try
            new_project(study;
                        data_dir     = ServerState.data_dir(),
                        projects_dir = ServerState.projects_dir(),
                        config_dir   = joinpath(dirname(ServerState.data_dir()), "config"))
            push!(initialised, study)
        catch e
            @warn "server: could not initialise project '$study'" exception=e
        end
    end
    initialised
end

## routes

@post "/api/v1/init" function(req)
    initialised = _ensure_all_projects()
    json((; initialised))
end

@get "/api/v1/studies" function(req)
    studies = map(_study_names()) do name
        runs   = _run_names(name)
        groups = _group_names(name)
        (; name,
           run_count        = length(runs),
           group_count      = length(groups),
           active_job_count = _active_job_count(name))
    end
    json(studies)
end

@get "/api/v1/studies/{study}" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    json(JSON3.read(_study_json(study)))
end
