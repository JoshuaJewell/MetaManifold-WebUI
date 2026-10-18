# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## R runtime
# RCall hosts exactly one embedded R session per Julia process, and R is neither
# thread-safe nor reentrant. Every caller that evaluates R must therefore
# serialise on the single lock defined here: the DADA2 pipeline stages, which run
# on spawned job threads, and the NMDS, PERMANOVA, and alpha-significance
# analyses, which run on HTTP handler tasks.
#
# One lock guards the interpreter. The pipeline resets the workspace with
# `rm(list=ls())` on acquire, which would destroy an analysis's globals
# mid-evaluation, and a pipeline runs on a background thread while the browser
# stays live.
module RRuntime

import RCall

export RBusyError, with_r_lock, r_busy, REQUIRED_R_PACKAGES, missing_r_packages

const _R_LOCK = ReentrantLock()

# The R packages the app loads, each with what stops working without it. install.jl
# checks the same list after restoring (test_r_runtime.jl keeps the two equal).
const REQUIRED_R_PACKAGES = [
    "dada2"         => "the DADA2 stages that run on this machine",
    "vegan"         => "NMDS, PERMANOVA and PERMDISP",
    "glmGamPoi"     => "shared dispersion-trend differential abundance",
    "MASS"          => "negative-binomial differential abundance",
    "zCompositions" => "CLR differential abundance (clr_welch)",
]

# Thrown when a caller that declines to wait indefinitely gives up on the runtime.
struct RBusyError <: Exception
    waited :: Float64
end

Base.showerror(io::IO, e::RBusyError) = print(io,
    "R runtime busy: another task, typically a pipeline run, held it for " *
    "longer than $(e.waited)s")

# True while any task holds the R runtime.
r_busy() = islocked(_R_LOCK)

"""
    with_r_lock(f; timeout=nothing) -> f()

Run `f` with exclusive access to the embedded R runtime.

With no `timeout` the caller waits for as long as it takes. That is what a
pipeline stage wants: it is the work the user explicitly asked for, and dropping
it because an analysis happened to be running would be worse than waiting.

With a `timeout` in seconds the caller waits at most that long and then throws
`RBusyError`. Interactive request handlers pass one, so that a pipeline run,
which can hold the runtime for hours, cannot hang an HTTP response indefinitely.
"""
function with_r_lock(f::Function; timeout::Union{Real,Nothing}=nothing)
    isnothing(timeout) && return lock(f, _R_LOCK)
    deadline = time() + timeout
    while true
        if trylock(_R_LOCK)
            try
                return f()
            finally
                unlock(_R_LOCK)
            end
        end
        time() >= deadline && throw(RBusyError(Float64(timeout)))
        sleep(0.05)
    end
end

"""
    missing_r_packages(packages=first.(REQUIRED_R_PACKAGES); timeout=nothing) -> Vector{String}

The named packages that the embedded R cannot load, in the order given.

Each is loaded, not merely looked up: a package can be listed in the library yet
fail to load, as when renv's links into an emptied cache dangle or a build no
longer matches the installed R.
"""
function missing_r_packages(packages::AbstractVector{<:AbstractString} = first.(REQUIRED_R_PACKAGES);
                            timeout::Union{Real,Nothing} = nothing)
    isempty(packages) && return String[]
    with_r_lock(; timeout) do
        RCall.globalEnv[:mm_required_pkgs] = String.(packages)
        try
            missing = RCall.rcopy(RCall.reval("""
                paste(mm_required_pkgs[!vapply(mm_required_pkgs, function(p)
                    suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE))),
                    logical(1))], collapse = "\\n")"""))
            isempty(missing) ? String[] : String.(split(missing, '\n'))
        finally
            RCall.reval("rm(mm_required_pkgs)")
        end
    end
end

end # module RRuntime
