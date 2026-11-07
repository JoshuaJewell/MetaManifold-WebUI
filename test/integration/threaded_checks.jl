# SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
# SPDX-License-Identifier: MPL-2.0

# @test inside Threads.@threads, with results that reach the enclosing testset.
#
# On Julia 1.12 the testset stack is task-local and a @threads task starts with an
# empty one, so a @test there reports to no testset: a pass is not counted, and a
# failure throws out of the task, ending that task's checks and every check after the
# loop in the enclosing testset, which reports one Error for the lot. From 1.13 tasks
# inherit the testset, so the same loop counts correctly and only CI's totals are short.
#
# threaded_checks runs one body per item in parallel, as @threads does. Each body gets
# a testset of its own to report to, and the results are replayed from the calling
# task once all of them are done, under one testset per item. The same code therefore
# counts every check on either version, keeps one item's failure from cutting the
# others short, and prints failures whole rather than interleaved.

using Test

struct CollectedChecks <: Test.AbstractTestSet
    description::String
    results::Vector{Any}
    CollectedChecks(description::AbstractString; _...) = new(description, Any[])
end
Test.record(ts::CollectedChecks, result) = (push!(ts.results, result); result)
Test.finish(ts::CollectedChecks) = ts

"""
    threaded_checks(body, items; name = string)

Call `body(i, item)` for each item, in parallel, reporting each item's `@test`s to the
enclosing testset under a nested testset called `name(item)`. An exception out of
`body` is recorded as that item's error and does not stop the others.
"""
function threaded_checks(body, items; name = string)
    collected = Vector{CollectedChecks}(undef, length(items))
    Threads.@threads for i in eachindex(items)
        collected[i] = @testset CollectedChecks "$(name(items[i]))" begin
            body(i, items[i])
        end
    end
    for checks in collected
        @testset "$(checks.description)" begin
            foreach(r -> Test.record(Test.get_testset(), r), checks.results)
        end
    end
    nothing
end
