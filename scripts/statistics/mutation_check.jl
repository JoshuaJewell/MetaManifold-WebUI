#!/usr/bin/env julia
# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Negative controls: plant a defect, prove a named test kills it.
#
# A passing test suite is evidence that the code currently agrees with its
# assertions. It is not evidence that the assertions can fail. This script makes
# the difference checkable by mutating one reviewed decision at a time and
# requiring the focused test file to go red.
#
# Every mutant is applied to the real tree, the focused tests are run, and the
# file is restored from the saved copy in a finally block. The script refuses to
# start when a target file is already modified, so a failed run cannot launder
# someone else's working-tree change.
#
# Usage:
#   julia --project=. scripts/statistics/mutation_check.jl \
#       [--julia julia] [--max-mutants N] [--output FILE] [--keep-logs DIR]
#
# Exit status:
#   0  every planted mutant was killed and the tree was restored
#   1  a mutant survived (the tests passed while the defect was present)
#   2  the campaign could not run (missing tool, dirty target, bad argument)

using SHA
using JSON3
using Dates
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

struct Mutant
    name::String
    path::String
    old::String
    new::String
    tests::Vector{String}
    description::String
end

const MUTANTS = Mutant[
    Mutant("catalogue-status-draft",
        "config/defaults/validated_statistics.yml",
        "status: \"pending_review\"", "status: \"approved\"",
        ["method_catalogue"],
        "A fork approval must not impersonate upstream scientific/release approval."),
    Mutant("catalogue-cell-limit-drift",
        "config/defaults/validated_statistics.yml",
        "max_cells: 100000", "max_cells: 100001",
        ["method_catalogue"],
        "A catalogue that advertises more than the code enforces must fail."),
    Mutant("request-limit-drift",
        "src/analysis/numeric_policy.jl",
        "const MAX_REQUEST_BYTES = 16 * 1024 * 1024",
        "const MAX_REQUEST_BYTES = 32 * 1024 * 1024",
        ["method_catalogue"],
        "The request limit is owned by one constant; drifting it must fail."),
    Mutant("cancellation-cleared-on-read",
        "src/analysis/numeric_policy.jl",
        "    token.cancelled[] &&",
        "    Threads.atomic_xchg!(token.cancelled, false) &&",
        ["statistics_adversarial"],
        "A second cancellation check must also refuse work."),
    Mutant("exact-proportion-becomes-float",
        "src/analysis/exact_summaries.jl",
        "proportions[i] = exact_relative_abundance(count, total; policy)",
        "proportions[i] = exact_relative_abundance(BigInt(round(Int, 1_000_000 * count / total)), 1_000_000; policy)",
        ["exact_summaries"],
        "Replacing an exact proportion with a float-derived one must be caught."),
]

_git_status_command(paths::Vector{String}) = `git -C $ROOT status --porcelain -- $paths`

function _git_status(paths::Vector{String})
    output = IOBuffer()
    try
        run(pipeline(_git_status_command(paths);
                     stdout=output, stderr=devnull))
    catch
        return "git status failed"
    end
    strip(String(take!(output)))
end

function _sha256(path::AbstractString)
    bytes2hex(sha256(read(path)))
end

_test_command(julia::String, tests::Vector{String}) =
    Cmd(`$julia --project=. -t2 --compiled-modules=existing test/runtests.jl $tests`; dir=ROOT)

function _run_tests(julia::String, tests::Vector{String}, log::String)
    command = _test_command(julia, tests)
    open(log, "w") do io
        println(io, "\$ ", command)
        flush(io)
        return run(pipeline(ignorestatus(command); stdout=io, stderr=io)).exitcode
    end
end

# Invalid catalogue mutations are rejected on loading, before an assertion can
# run. Name those exact expected errors; an arbitrary exception is not a kill.
const EXPECTED_REJECTIONS = Dict(
    "catalogue-status-draft" => "catalogue must retain pending_review until independent review and release approval",
    "catalogue-cell-limit-drift" => "catalogue.limits.max_cells is 100001 but the code enforces 100000",
    "request-limit-drift" => "catalogue.limits.max_request_bytes is 16777216 but the code enforces 33554432",
)

# A crash, parse error, missing dependency or killed process is not a detected
# defect. Require a completed test report and an assertion or named rejection.
function _expected_test_failure(exitcode::Integer, transcript::AbstractString;
                                expected_rejection::Union{Nothing,String}=nothing)
    exitcode == 1 && occursin("Test Summary:", transcript) || return false
    if occursin("Error During Test at ", transcript)
        return !isnothing(expected_rejection) &&
            occursin("ArgumentError: " * expected_rejection, transcript)
    end
    occursin("Test Failed at ", transcript)
end

function _argument_value(args, index, name)
    index < length(args) || error("missing value for $name")
    args[index + 1]
end

function main(args::Vector{String})
    julia = "julia"
    max_mutants = length(MUTANTS)
    output = "mutation-report.json"
    keep_logs = nothing
    index = 1
    while index <= length(args)
        argument = args[index]
        if argument == "--julia"
            julia = _argument_value(args, index, argument); index += 2
        elseif argument == "--max-mutants"
            max_mutants = parse(Int, _argument_value(args, index, argument)); index += 2
        elseif argument == "--output"
            output = _argument_value(args, index, argument); index += 2
        elseif argument == "--keep-logs"
            keep_logs = _argument_value(args, index, argument); index += 2
        else
            println(stderr, "mutation_check: unknown argument '$argument'")
            return 2
        end
    end
    1 <= max_mutants <= length(MUTANTS) || (println(stderr,
        "mutation_check: --max-mutants must be in 1:$(length(MUTANTS))"); return 2)

    selected = MUTANTS[1:max_mutants]
    paths = sort(unique([mutant.path for mutant in selected]))
    dirty = _git_status(paths)
    isempty(dirty) || (println(stderr,
        "mutation_check: refusing to run on a modified target:\n$dirty"); return 2)
    # Check every anchor before spending time on the baseline or mutating files.
    for mutant in selected
        original = read(joinpath(ROOT, mutant.path), String)
        occurrences = length(split(original, mutant.old)) - 1
        occurrences == 1 || (println(stderr,
            "mutation_check: '$(mutant.name)' anchor matched $occurrences times; refusing"); return 2)
    end

    logs = isnothing(keep_logs) ? mktempdir(; cleanup=false) : String(keep_logs)
    mkpath(logs)
    baseline_tests = sort(unique(vcat([mutant.tests for mutant in selected]..., ["statistics_evidence"])))
    baseline_log = joinpath(logs, "00-pristine.log")
    _run_tests(julia, baseline_tests, baseline_log) == 0 || (println(stderr,
        "mutation_check: pristine tests failed; no mutant can be counted as killed; see $baseline_log"); return 2)
    results = Vector{Dict{String,Any}}()
    survived = String[]
    for (index, mutant) in enumerate(selected)
        target = joinpath(ROOT, mutant.path)
        original = read(target, String)
        occurrences = length(split(original, mutant.old)) - 1
        occurrences == 1 || (println(stderr,
            "mutation_check: '$(mutant.name)' anchor matched $occurrences times " *
            "in $(mutant.path); refusing"); return 2)
        backup = joinpath(logs, @sprintf("%02d-%s.orig", index, mutant.name))
        log = joinpath(logs, @sprintf("%02d-%s.log", index, mutant.name))
        write(backup, original)
        killed = false
        exitcode = -1
        try
            write(target, replace(original, mutant.old => mutant.new))
            exitcode = _run_tests(julia, mutant.tests, log)
            killed = _expected_test_failure(exitcode, read(log, String);
                expected_rejection=get(EXPECTED_REJECTIONS, mutant.name, nothing))
        finally
            write(target, original)
        end
        _sha256(target) == _sha256(backup) || (println(stderr,
            "mutation_check: failed to restore $(mutant.path)"); return 2)
        exitcode != 0 && !killed && (println(stderr,
            "mutation_check: invalid test failure for $(mutant.name); see $log"); return 2)
        killed || push!(survived, mutant.name)
        push!(results, Dict{String,Any}(
            "name" => mutant.name,
            "description" => mutant.description,
            "path" => mutant.path,
            "tests" => mutant.tests,
            "killed" => killed,
            "test_exitcode" => exitcode,
            "expected_rejection" => get(EXPECTED_REJECTIONS, mutant.name, nothing),
            "log" => log,
        ))
        println(killed ? "KILLED   " : "SURVIVED ", mutant.name)
    end

    isempty(_git_status(paths)) || (println(stderr,
        "mutation_check: the tree was not restored"); return 2)

    report = Dict{String,Any}(
        "schema" => "metamanifold.statistics-mutation-report.v1",
        "command" => "scripts/statistics/mutation_check.jl",
        "julia" => string(VERSION),
        "generated" => string(Dates.now(UTC)),
        "mutants" => results,
        "killed" => count(result -> result["killed"], results),
        "survived" => survived,
        "complete" => isempty(survived),
    )
    open(output, "w") do io
        JSON3.pretty(io, report)
        println(io)
    end
    println("mutation_check: $(report["killed"])/$(length(results)) mutants killed; ",
            "report $output")
    return isempty(survived) ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
