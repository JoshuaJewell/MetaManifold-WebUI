# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

using Test
using MetaManifold.NumericPolicy
using MetaManifold.HighPrecision
using MetaManifold.Canonical

# Independently produced with Python decimal. These are reference values, not
# values copied from Julia/MPFR; REF_LN2 has enough digits for the 1024-bit test.
const HP_FIXTURE = [1 2; 1 0; 1 0]
const HP_LABELS = ["A", "B"]
const HP_FEATURES = ["f1", "f2", "f3"]

const REF_LN2 = "0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202357581305570326707516350759619307275708283714351903070386238916734711233501153644979552391204751726815749320651555247341395258829504530070953263666426541042391578149520437404303855008019441706416715186447128399681717845469570262716310645461502572074024816377733896385506952606683411372738737229289564935470257626520988596932019650586"
const REF_LN3 = "1.098612288668109691395245236922525704647490557822749451734694333637494293218608966873615754813732088787970029065957865742368004225930519821052801870767277410603162769183381367179373698844360959903742570316795911521145591917750671347054940166775580222203170252946897560690106521505642868138036317373298577782366991654792131818149020030103823630122248652748198225991097452490896458053467008845965085748444119018857087647494867079613085829411602166121184001409826"
const REF_TWO_THIRDS = "0.6666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667"
const REF_H_1234 = "1.279854225833667467182185762149429928131354186921692886392239591382930699742075074568049151649059984163575294301928203177013603961545520253317776056233522504371230806742168835499680746549033237086147384098562106321046048112868339308649672681189135461287570216937288791515535715833594380020178552966375097993963900928200120047744717834903206362014300458629735272704036751340598445797090470318251133760850603018005627582297102594277585911107434593371009711254916"

function test_bigfloat_close(actual, reference, bits, ulps=16)
    setprecision(BigFloat, bits) do
        expected = BigFloat(reference; precision=bits)
        one = BigFloat(1; precision=bits)
        tolerance = max(abs(expected), one) * ulps /
                    (BigFloat(2; precision=bits)^(bits - 4))
        @test isfinite(actual)
        @test abs(actual - expected) <= tolerance
    end
end

@testset "High-precision diversity known answers" begin
    policy = numeric_policy(mode=:high_precision, precision_bits=512, round_digits=30)
    s = high_precision_summary(HP_FIXTURE; sample_labels=HP_LABELS,
        feature_labels=HP_FEATURES, policy=policy)
    @test s.achieved_bits == 512
    @test s.samples[1].richness == 3
    @test s.samples[2].richness == 1
    @test s.samples[1].total == 3
    test_bigfloat_close(s.samples[1].shannon, REF_LN3, 512)
    test_bigfloat_close(s.samples[1].gini_simpson, REF_TWO_THIRDS, 512)
    test_bigfloat_close(s.samples[2].shannon, "0", 512)
    test_bigfloat_close(s.samples[2].gini_simpson, "0", 512)

    other = high_precision_summary(reshape([1, 2, 3, 4], 4, 1);
        sample_labels=["counts"], feature_labels=["a", "b", "c", "d"], policy=policy)
    @test other.samples[1].richness == 4
    test_bigfloat_close(other.samples[1].shannon, REF_H_1234, 512)
    test_bigfloat_close(other.samples[1].gini_simpson, "0.7", 512)

    storage = high_precision_to_storage(other)
    @test storage["status"] == "success"
    @test storage["input_status"] == "exact"
    @test storage["achieved_precision_bits"] == 512
    @test storage["samples"][1]["richness"] == 4
    @test storage["samples"][1]["shannon"]["kind"] == "approximate"
    @test storage["samples"][1]["shannon"]["precision_bits"] == 512
    @test storage["samples"][1]["shannon"]["display_kind"] == "rounded"
    @test length(jcs_sha256(storage)) == 64
    @test occursin("not exact", high_precision_to_display(other))

    encoded = storage["samples"][1]["shannon"]["value"]
    round_trip = parse_approximate_decimal(encoded; bits=512)
    @test round_trip.value == other.samples[1].shannon
end

@testset "High precision does not upgrade floats or define empty samples" begin
    policy = numeric_policy(mode=:high_precision, precision_bits=256)
    @test_throws UnsupportedRepresentationError high_precision_summary(Float64.(HP_FIXTURE);
        sample_labels=HP_LABELS, feature_labels=HP_FEATURES, policy=policy)
    @test_throws ArgumentError high_precision_summary([1 -2; 3 4];
        sample_labels=HP_LABELS, feature_labels=["x", "y"], policy=policy)
    # A policy for another mode is a caller error, not a value that cannot be
    # represented: assert_mode raises ArgumentError.
    @test_throws ArgumentError high_precision_summary(HP_FIXTURE;
        sample_labels=HP_LABELS, feature_labels=HP_FEATURES,
        policy=numeric_policy(:exact_counts))

    empty = high_precision_summary(zeros(Int, 3, 1); sample_labels=["empty"],
        feature_labels=HP_FEATURES, policy=policy)
    @test empty.samples[1].richness == 0
    @test isnothing(empty.samples[1].shannon)
    @test isnothing(empty.samples[1].gini_simpson)
    storage = high_precision_to_storage(empty)
    @test storage["status"] == "partial"
    @test storage["samples"][1]["shannon"] === nothing
    @test storage["samples"][1]["status"] == "undefined_zero_total"
    @test occursin("undefined (zero total)", high_precision_to_display(empty))
end

@testset "Extreme scale retains nonzero high-precision information" begin
    # A Float64 conversion underflows the rare-taxon Shannon contribution to 0;
    # MPFR retains it. The rational budget must accommodate N^2 for Gini-Simpson.
    counts = reshape(BigInt[big(2)^4_000, 1], 2, 1)
    policy = numeric_policy(mode=:high_precision, precision_bits=256,
                            max_denominator_bits=50_000, round_digits=12)
    s = high_precision_summary(counts; sample_labels=["astronomical"],
        feature_labels=["common", "rare"], policy=policy)
    @test s.samples[1].richness == 2
    @test s.samples[1].shannon > 0
    @test Float64(s.samples[1].shannon) == 0.0
    @test 0 < s.samples[1].gini_simpson < 1
    storage = high_precision_to_storage(s)
    @test storage["samples"][1]["total"] == string(big(2)^4_000 + 1)
    @test parse_approximate_decimal(storage["samples"][1]["shannon"]["value"]; bits=256).value ==
          s.samples[1].shannon
end

@testset "Concurrent precision scopes do not leak" begin
    old = precision(BigFloat)
    try
        setprecision(BigFloat, 192)
        levels = [128, 256, 512, 1024, 256, 512, 128, 1024]
        results = Vector{Any}(undef, length(levels))
        @sync for i in eachindex(levels)
            Threads.@spawn begin
                p = numeric_policy(mode=:high_precision, precision_bits=levels[i])
                results[i] = high_precision_summary(reshape([1, 1], 2, 1);
                    sample_labels=["s$i"], feature_labels=["a", "b"], policy=p)
            end
        end
        @test precision(BigFloat) == 192
        for (i, result) in enumerate(results)
            @test result.achieved_bits == levels[i]
            @test precision(result.samples[1].shannon) == levels[i]
            @test precision(result.samples[1].gini_simpson) == levels[i]
            test_bigfloat_close(result.samples[1].shannon, REF_LN2, levels[i], 32)
        end
    finally
        setprecision(BigFloat, old)
    end
end
