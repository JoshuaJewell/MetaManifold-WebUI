# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/jobs
using JSON3

function _job_to_namedtuple(j::Job)
    (; j.id, j.type,
       j.study, j.run, j.group,
       j.stage, j.message,
       status      = string(j.status),
       created_at  = string(j.created_at),
       finished_at = isnothing(j.finished_at) ? nothing : string(j.finished_at))
end

@get "/api/v1/jobs" function(req)
    study  = get(queryparams(req), "study",  nothing)
    status = get(queryparams(req), "status", nothing)
    jobs   = list_jobs(; study, status)
    json(map(_job_to_namedtuple, jobs))
end

@get "/api/v1/jobs/{id}" function(req, id::String)
    j = get_job(id)
    isnothing(j) && return json_error(404, "job_not_found", "Job '$id' not found")
    json(_job_to_namedtuple(j))
end

@delete "/api/v1/jobs/{id}" function(req, id::String)
    # Deleting a finished job succeeds.
    isnothing(get_job(id)) && return json_error(404, "job_not_found", "Job '$id' not found")
    cancel_job!(id)
    HTTP.Response(204)
end
