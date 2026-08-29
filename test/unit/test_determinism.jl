# Determinism of every output this codebase computes from a results table, from
# config, and from the stochastic normalisation / ordination steps.
#
# The central device is an adversarial engine matrix. The same logical table is
# built under every combination of
#
#   DuckDB threads              the process's thread count, and that plus 12
#   preserve_insertion_order    true, false
#   physical row order          sorted, shuffled
#
# and each output is required to be byte-identical across all eight. A guard test
# first shows the matrix really does perturb unordered queries, so the invariance
# checks cannot pass vacuously. The table is large enough (> 2 DuckDB row groups)
# that scans and aggregations actually parallelise.
#
# What this file cannot reach - DADA2's own thread/seed behaviour, vsearch and
# cd-hit threading, and the local-versus-remote R stacks - is exercised by
# test/determinism/ on the compute server.

using MetaManifold
SV = MetaManifold.Server
_saved_root = SV.ServerState._root[]

using Random: MersenneTwister, randperm, shuffle
using SHA: sha256
using OrderedCollections: OrderedDict

const _DET_N       = 300_000                     # ~2.5 DuckDB row groups
const _DET_SAMPLES = ["S01", "S02", "S03", "S04", "S05", "S06"]
const _DET_GENERA  = ["Blastocystis", "Giardia", "Entamoeba", "Hexamastix",
                      "Escherichia", "Clostridium", "Candida", ""]

# The table is generated inside DuckDB from `range()`, with every value a pure
# function of the row index, so each configuration holds the same logical rows.
# Building it with one CREATE TABLE AS keeps the suite fast; appending 300k rows
# from Julia one value at a time, twenty times over, did not.
#
# Two rows carry engineered ties on the aggregated read total: genera TieAlpha
# and TieBeta get identical per-sample sums, so only a tiebreak orders them.
function _det_db(; threads::Int, preserve::Bool, shuffled::Bool)
    db  = DuckDB.DB()
    con = DBInterface.connect(db)
    DBInterface.execute(con, "SET threads = $threads")
    # Keep insertion order while building so the chosen physical order sticks;
    # the setting under test is applied afterwards, for the queries.
    DBInterface.execute(con, "SET preserve_insertion_order = true")
    genera = join(["'$g'" for g in _DET_GENERA], ", ")
    count_col(k) = """CASE WHEN i IN (1, 2) THEN $(k == 1 ? 1000 : 0)
                           WHEN hash(i * 101 + $k) % 100 < 35
                           THEN CAST(1 + hash(i * 211 + $k) % 40 AS BIGINT)
                           ELSE 0 END"""
    order = shuffled ? "ORDER BY hash(i * 7919 + 13)" : "ORDER BY i"
    DBInterface.execute(con, """
        CREATE TABLE merged AS
        SELECT 'seq' || i                          AS "SeqName",
               'ACGT' || i                         AS "sequence",
               CAST(97.0 AS DOUBLE)                AS "Pident",
               CASE WHEN i = 1 THEN 'TieAlpha' WHEN i = 2 THEN 'TieBeta'
                    ELSE [$genera][((i - 1) % $(length(_DET_GENERA))) + 1] END AS "Genus",
               CAST(hash(i * 307) % 101 AS DOUBLE) AS "Genus_boot",
               $(join(["$(count_col(k)) AS \"$(_DET_SAMPLES[k])\"" for k in 1:6], ",\n               "))
        FROM range(1, $(_DET_N + 1)) t(i)
        $order""")
    DBInterface.execute(con, "SET preserve_insertion_order = $preserve")
    db, con
end

# DuckDB.jl registers every Julia thread as an external DuckDB thread and refuses
# a thread count below that, so the low end is the process's own thread count.
const _DET_THREADS = (Threads.nthreads(), Threads.nthreads() + 12)
const _DET_MATRIX = [(; threads, preserve, shuffled)
                     for threads in _DET_THREADS, preserve in (true, false), shuffled in (false, true)]

# Fingerprints hash raw bytes, not `repr` text: the matrices here run to
# millions of cells, and rendering them as strings eight times over dominated
# the suite's runtime. Shape is folded in so a reshaped array cannot collide.
# Hash a byte vector, never the String itself: SHA.jl's update! checks aliasing
# with `objectid`, which for a String re-hashes the entire contents on every
# 64-byte block, making a multi-megabyte fingerprint quadratic.
_det_hash(x::AbstractString) = bytes2hex(sha256(Vector{UInt8}(codeunits(x))))
_det_hash(x::AbstractArray{<:Real}) =
    bytes2hex(sha256(vcat(Vector{UInt8}(string(size(x))), reinterpret(UInt8, vec(Float64.(x))))))
_det_hash(x::AbstractVector{<:AbstractString}) = _det_hash(join(x, '\x1f'))
_det_hash(x::Tuple) = _det_hash(join(map(_det_hash, x), ","))
function _det_hash(df::DataFrame)
    io = IOBuffer(); CSV.write(io, df)
    _det_hash(String(take!(io)))
end
_det_hash(x) = _det_hash(repr(x))
_det_json(resp) = String(resp.body)
_det_body(d) = JSON3.read(JSON3.write(d))

# Every output under test, computed against one connection. Returns an ordered
# name => fingerprint map so a failure names the output that diverged.
function _det_outputs(con)
    out = OrderedDict{String,String}()
    table   = "merged"
    columns = SV._duckdb_columns(con, table)
    scols   = SV.Analysis.sample_columns(con, table)
    out["sample_columns"] = _det_hash(scols)

    params = SV._body_filter_params(_det_body(Dict(
        "colFilters" => Dict("Genus_boot" => Dict("min" => 20),
                             SV.SAMPLE_READS_FILTER_KEY => Dict("min" => 1)))))
    where, wp = SV._analysis_where_clause(params, columns)
    kept = SV._retain_sample_columns(con, table, scols, params, where, wp)
    out["retained_samples"] = _det_hash(kept)

    mat = SV.Analysis.filtered_counts(con, table, kept, where, wp)
    out["filtered_counts"] = _det_hash(mat)

    for method in ("rarefy", "srs")
        norm = SV.normalise_counts(mat; method, depth=0, seed=123)
        out["normalise.$method"] = _det_hash((norm.mat, norm.kept))
        alpha = [(SV.richness(round.(Int, norm.mat[i, :])),
                  SV.shannon(round.(Int, norm.mat[i, :])),
                  SV.simpson(round.(Int, norm.mat[i, :]))) for i in axes(norm.mat, 1)]
        out["alpha.$method"] = _det_hash(alpha)
    end

    agg = SV.Analysis.aggregate_by_taxon(con, table, kept, "Genus", where, wp)
    out["aggregate_by_taxon"] = _det_hash(agg)

    out["venn_taxa_present"] = _det_hash(SV.Analysis.venn_taxa_present(con, table, kept, "Genus", where, wp))

    cd = SV._chart_data(con, table, columns, kept, where, wp, "rank", "Genus",
                        nothing, String[], "run", 5, "VSEARCH")
    for relative in (true, false)
        fig = SV.bar_chart(cd.segment_labels, cd.sample_names, cd.counts;
                           top_n=cd.effective_top_n, relative, mode="stacked")
        out["bar_chart.relative=$relative"] = _det_hash(JSON3.write(fig))
    end

    # Pagination: unsorted, sorted on a heavily tied column, and descending.
    for (label, extra) in (("unsorted", Dict()),
                           ("genus_asc", Dict("sortBy" => "Genus")),
                           ("boot_desc", Dict("sortBy" => "Genus_boot", "sortDir" => "desc")))
        pages = String[]
        for page in (1, 2, 37)
            body = _det_body(merge(Dict("page" => page, "perPage" => 997,
                "colFilters" => Dict(SV.SAMPLE_READS_FILTER_KEY => Dict("min" => 1)),
                "filter" => ""), extra))
            push!(pages, _det_json(SV._duckdb_paginated_query(con, table, body)))
        end
        out["paginated.$label"] = _det_hash(join(pages, "\n"))
    end

    out["save_export_frame"] = _det_hash(SV._filtered_table_df(con, table, columns, params, "Genus", "asc"))

    out["filtered_df"] = _det_hash(SV.Analysis.filtered_df(con, table, where, wp))
    out
end

@testset "Determinism" begin

    @testset "the engine matrix really perturbs unordered queries (non-vacuity guard)" begin
        firsts = Set{String}()
        group_orders = Set{String}()
        for cfg in _DET_MATRIX
            db, con = _det_db(; cfg...)
            try
                r = DataFrame(DBInterface.execute(con, "SELECT \"SeqName\" FROM merged LIMIT 50"))
                push!(firsts, join(r.SeqName, ","))
                g = DataFrame(DBInterface.execute(con,
                    "SELECT \"Genus\", SUM(\"S01\") AS n FROM merged GROUP BY \"Genus\""))
                push!(group_orders, join(g.Genus, ","))
            finally
                DBInterface.close!(con); close(db)
            end
        end
        # Unordered scans differ with physical order; unordered GROUP BY output
        # order is not fixed either. If neither varied, the matrix would prove
        # nothing about the ORDER BYs under test.
        @test length(firsts) > 1
        @info "Engine matrix: $(length(firsts)) distinct unordered scan orders, " *
              "$(length(group_orders)) distinct unordered GROUP BY orders across $(length(_DET_MATRIX)) configurations"
    end

    @testset "table-derived outputs are identical across the engine matrix" begin
        baseline = nothing
        base_cfg = nothing
        for cfg in _DET_MATRIX
            db, con = _det_db(; cfg...)
            outs = try
                _det_outputs(con)
            finally
                DBInterface.close!(con); close(db)
            end
            if isnothing(baseline)
                baseline, base_cfg = outs, cfg
                continue
            end
            for (name, h) in outs
                ok = h == baseline[name]
                ok || @error "Output diverged across engine configurations" output=name baseline=base_cfg this=cfg
                @test ok
            end
        end
        # Cross-process check: MM_DET_DUMP=<path> writes the fingerprints so two
        # separate processes (different Julia thread counts, a fresh session) can
        # be compared byte for byte - see test/determinism/cross_process.sh.
        dump = get(ENV, "MM_DET_DUMP", "")
        if !isempty(dump)
            open(dump, "w") do io
                for (name, h) in baseline
                    println(io, name, "\t", h)
                end
            end
        end
        # And the same outputs repeated on one connection.
        db, con = _det_db(threads=last(_DET_THREADS), preserve=false, shuffled=true)
        try
            @test _det_outputs(con) == _det_outputs(con)
        finally
            DBInterface.close!(con); close(db)
        end
    end

    @testset "sample_columns follows the table's column order" begin
        db = DuckDB.DB(); con = DBInterface.connect(db)
        try
            # Deliberately non-alphabetical, interleaved with non-count columns.
            DBInterface.execute(con, """CREATE TABLE t ("Zeta" BIGINT, "SeqName" VARCHAR,
                "alpha" BIGINT, "Genus_boot" DOUBLE, "Mid" DOUBLE, "Genus" VARCHAR, "b2" INTEGER)""")
            @test SV.Analysis.sample_columns(con, "t") == ["Zeta", "alpha", "Mid", "b2"]
        finally
            DBInterface.close!(con); close(db)
        end
    end

    @testset "aggregate_by_taxon breaks read-total ties by name" begin
        db, con = _det_db(threads=last(_DET_THREADS), preserve=false, shuffled=true)
        try
            agg = SV.Analysis.aggregate_by_taxon(con, "merged", ["S01"], "Genus",
                "WHERE \"Genus\" IN ('TieAlpha','TieBeta')", Any[])
            @test agg.taxon == ["TieAlpha", "TieBeta"]
        finally
            DBInterface.close!(con); close(db)
        end
    end

    @testset "pagination is a partition: every row exactly once, even on tied sorts" begin
        db, con = _det_db(threads=last(_DET_THREADS), preserve=false, shuffled=true)
        try
            DBInterface.execute(con, "CREATE TABLE small AS SELECT * FROM merged WHERE \"SeqName\" LIKE 'seq1%' AND length(\"SeqName\") <= 5")
            n = only(DataFrame(DBInterface.execute(con, "SELECT COUNT(*) AS n FROM small"))).n
            for sort_body in (Dict(), Dict("sortBy" => "Genus"), Dict("sortBy" => "Genus_boot", "sortDir" => "desc"))
                seen = String[]
                per = 37
                for page in 1:cld(n, per)
                    d = JSON3.read(_det_json(SV._duckdb_paginated_query(con, "small",
                        _det_body(merge(Dict("page" => page, "perPage" => per), sort_body)))))
                    append!(seen, [string(r.SeqName) for r in d.rows])
                end
                @test length(seen) == n
                @test allunique(seen)
            end
        finally
            DBInterface.close!(con); close(db)
        end
    end

    @testset "normalisation: seeded, and sensitive to what it must be sensitive to" begin
        rng = MersenneTwister(1)
        mat = Float64.(rand(rng, 0:30, 6, 400))
        for method in ("rarefy", "srs")
            a = SV.normalise_counts(mat; method, depth=0, seed=42)
            b = SV.normalise_counts(copy(mat); method, depth=0, seed=42)
            @test a.mat == b.mat && a.kept == b.kept
            # Inputs are not mutated (a mutation would make a second call differ).
            @test mat == Float64.(rand(MersenneTwister(1), 0:30, 6, 400))
        end
        # Rarefaction genuinely depends on the seed, so a fixed seed is load-bearing.
        @test SV.rarefy(mat; depth=100, seed=1) != SV.rarefy(mat; depth=100, seed=2)
        # ...and on feature order: the same data with its columns permuted and
        # un-permuted afterwards gives different draws. This is why the count
        # matrix's feature order has to be pinned by the query.
        p = randperm(MersenneTwister(3), size(mat, 2))
        diffs = count(s -> SV.rarefy(mat[:, p]; depth=100, seed=s)[:, invperm(p)] !=
                           SV.rarefy(mat; depth=100, seed=s), 1:5)
        @test diffs > 0
        # ...and on sample (row) order, since one stream is consumed row by row.
        q = [2, 1, 3, 4, 5, 6]
        @test SV.rarefy(mat[q, :]; depth=100, seed=9)[q, :] != SV.rarefy(mat; depth=100, seed=9)
    end

    @testset "combined cross-run matrices ignore input row order" begin
        mk(rows) = DataFrame(sequence=[r[1] for r in rows], a=[r[2] for r in rows], b=[r[3] for r in rows])
        rows = [("AAA", 1.0, 2.0), ("CCC", 3.0, 0.0), ("GGG", 0.0, 5.0), ("TTT", 7.0, 7.0)]
        r1 = SV.combined_asv_counts_across_runs([("r1", ["a", "b"], mk(rows)), ("r2", ["a", "b"], mk(reverse(rows)))])
        r2 = SV.combined_asv_counts_across_runs([("r1", ["a", "b"], mk(shuffle(MersenneTwister(5), rows))),
                                                ("r2", ["a", "b"], mk(rows))])
        @test r1 == r2
        @test r1[3] == ["AAA", "CCC", "GGG", "TTT"]
    end

    @testset "YAML documents are written in canonical key order" begin
        keys_ = ["zeta", "alpha", "mid", "beta", "omega", "gamma"]
        d1 = Dict{String,Any}(); for k in keys_; d1[k] = Dict("y" => 1, "x" => [Dict("b" => 1, "a" => 2)]); end
        d2 = Dict{String,Any}(); sizehint!(d2, 1024)
        for k in reverse(keys_); d2[k] = Dict("x" => [Dict("a" => 2, "b" => 1)], "y" => 1); end
        # Grow and shrink to give d2 a different internal layout from d1.
        for i in 1:500; d2["tmp$i"] = i; end
        for i in 1:500; delete!(d2, "tmp$i"); end
        @test YAML.write(SV.Config.canonical_yaml_doc(d1)) == YAML.write(SV.Config.canonical_yaml_doc(d2))
        @test startswith(YAML.write(SV.Config.canonical_yaml_doc(d1)), "alpha:")
        # An OrderedDict's chosen order survives; sequences keep their order.
        od = OrderedDict("second" => 1, "first" => 2)
        @test startswith(YAML.write(SV.Config.canonical_yaml_doc(od)), "second:")
        @test SV.Config.canonical_yaml_doc(Any[3, 1, 2]) == Any[3, 1, 2]
    end

    @testset "run_config.yml is identical whatever the source files' key order" begin
        texts = String[]
        for variant in 1:2
            dir = mktempdir()
            try
                config_dir = joinpath(dir, "config")
                mkpath(joinpath(config_dir, "defaults"))
                default_body = variant == 1 ?
                    "seed: 123\ncutadapt:\n  min_length: 200\n  cores: 0\nanalysis:\n  normalisation: none\n  transform: none\n" :
                    "analysis:\n  transform: none\n  normalisation: none\ncutadapt:\n  cores: 0\n  min_length: 200\nseed: 123\n"
                write(joinpath(config_dir, "defaults", "pipeline.yml"), default_body)
                data_dir = joinpath(dir, "data", "study", "run")
                mkpath(data_dir)
                write(joinpath(dir, "data", "study", "pipeline.yml"),
                      variant == 1 ? "cutadapt:\n  min_length: 150\nswarm:\n  d: 1\n" :
                                     "swarm:\n  d: 1\ncutadapt:\n  min_length: 150\n")
                proj_dir = joinpath(dir, "projects", "study", "run")
                mkpath(proj_dir)
                ctx = SV.ProjectCtx(proj_dir, config_dir, data_dir,
                                    joinpath(dir, "projects", "study"), joinpath(dir, "data", "study"))
                push!(texts, read(SV.Config.write_run_config(ctx), String))
            finally
                rm(dir; recursive=true, force=true)
            end
        end
        @test texts[1] == texts[2]
    end

    @testset "attestation run/config sections are canonical" begin
        texts = String[]
        for (run, cfg) in ((Dict("study" => "s", "group" => nothing, "run" => "r"),
                            Dict{String,Any}("seed" => 1, "analysis" => Dict("b" => 1, "a" => 2))),
                           (Dict("run" => "r", "study" => "s", "group" => nothing),
                            Dict{String,Any}("analysis" => Dict("a" => 2, "b" => 1), "seed" => 1)))
            path = tempname() * ".yml"
            att = SV.MetaManifold.Provenance.Attestation(; run, config=cfg, config_sha256="00")
            SV.MetaManifold.Provenance.write_attestation(att, path; merge_existing=false)
            # The timestamp is the one field that is meant to change.
            push!(texts, join(filter(l -> !occursin("generated", l), readlines(path)), "\n"))
            rm(path; force=true)
        end
        @test texts[1] == texts[2]
    end

    @testset "vsearch hits are written in query order, keeping per-query ranking" begin
        td = mktempdir()
        try
            hits = ["seq3\tRefA\t99.0", "seq3\tRefB\t97.5", "seq10\tRefC\t88.0",
                    "seq1\tRefD\t91.0", "otu2\tRefE\t95.0", "seq1\tRefF\t90.0"]
            outs = String[]
            for seed in 1:4
                path = joinpath(td, "taxonomy_$seed.tsv")
                # Queries in thread-finish order, but each query's hits stay in the
                # rank order vsearch emitted them in, as vsearch writes them.
                blocks = [["seq3\tRefA\t99.0", "seq3\tRefB\t97.5"], ["seq10\tRefC\t88.0"],
                          ["seq1\tRefD\t91.0", "seq1\tRefF\t90.0"], ["otu2\tRefE\t95.0"]]
                write(path, join(vcat(shuffle(MersenneTwister(seed), blocks)...), "\n") * "\n")
                SV.Tools._sort_hits_by_query!(path)
                push!(outs, read(path, String))
            end
            @test allequal(outs)
            @test split(strip(outs[1]), "\n") ==
                  ["otu2\tRefE\t95.0", "seq1\tRefD\t91.0", "seq1\tRefF\t90.0",
                   "seq10\tRefC\t88.0", "seq3\tRefA\t99.0", "seq3\tRefB\t97.5"]
            @test isempty(filter(f -> startswith(f, "jl_"), readdir(td)))   # no temp files left behind
        finally
            rm(td; recursive=true, force=true)
        end
    end

    @testset "run_config.yml is never observed half-written by a concurrent reader" begin
        # Needs a reader running while the writer writes; with one thread the
        # reader task only runs between writes, so the check proves nothing.
        if Threads.nthreads() < 2
            @info "single-threaded process - skipping the run_config write race check"
            @test true
        else
        dir = mktempdir()
        try
            config_dir = joinpath(dir, "config")
            mkpath(joinpath(config_dir, "defaults"))
            # A config large enough that a non-atomic write is a window, not an instant.
            body = join(["section_$i:\n" * join(["  key_$j: $(i * j)" for j in 1:40], "\n") for i in 1:60], "\n") * "\n"
            write(joinpath(config_dir, "defaults", "pipeline.yml"), body)
            data_dir = joinpath(dir, "data", "study", "run"); mkpath(data_dir)
            study_yml = joinpath(dir, "data", "study", "pipeline.yml")
            proj_dir = joinpath(dir, "projects", "study", "run"); mkpath(proj_dir)
            ctx = SV.ProjectCtx(proj_dir, config_dir, data_dir,
                                joinpath(dir, "projects", "study"), joinpath(dir, "data", "study"))
            # The first version needs a seed too, or a reader that starts before
            # the loop's first rewrite sees a complete file without one.
            write(study_yml, "seed: 0\n")
            path = SV.Config.write_run_config(ctx)

            # Readers check raw bytes, not a YAML parse: a parse is slow enough that
            # readers barely overlap the writes and the test could not fail. Keys
            # are written sorted, so a complete file always ends with the seed line.
            complete(text) = occursin(r"\nseed: \d+\n\z", text)
            stop = Threads.Atomic{Bool}(false)
            bad  = Threads.Atomic{Int}(0)
            reads = Threads.Atomic{Int}(0)
            readers = [Threads.@spawn begin
                while !stop[]
                    text = try read(path, String) catch; "" end
                    complete(text) || Threads.atomic_add!(bad, 1)
                    Threads.atomic_add!(reads, 1)
                end
            end for _ in 1:max(1, Threads.nthreads() - 1)]
            for i in 1:400
                # Each edit makes the source newer, forcing a regeneration.
                write(study_yml, "seed: $i\n")
                SV.Config.write_run_config(ctx)
            end
            stop[] = true
            foreach(wait, readers)
            @test reads[] > 0
            @test bad[] == 0
            @test YAML.load_file(path)["seed"] == 400
            @test isempty(filter(f -> startswith(f, "jl_"), readdir(proj_dir)))
        finally
            rm(dir; recursive=true, force=true)
        end
        end
    end

    @testset "CD-HIT count collapse writes a stable, input-ordered table" begin
        td = mktempdir()
        try
            clstr = joinpath(td, "a.clstr")
            write(clstr, ">Cluster 0\n0\t300nt, >seq9... *\n1\t295nt, >seq2... at 98.00%\n" *
                         ">Cluster 1\n0\t250nt, >seq5... *\n>Cluster 2\n0\t250nt, >seq1... *\n" *
                         "1\t240nt, >seq7... at 97.00%\n")
            counts = joinpath(td, "c.csv")
            write(counts, "SeqName,Sequence,A,B\nseq9,AAA,1,2\nseq2,CCC,3,4\nseq5,GGG,5,6\n" *
                          "seq1,TTT,7,8\nseq7,ACG,9,10\n")
            outs = String[]
            for i in 1:3
                o = joinpath(td, "out$i.csv")
                SV.Tools._collapse_cdhit_counts(clstr, counts, o)
                push!(outs, read(o, String))
            end
            @test allequal(outs)
            df = CSV.read(joinpath(td, "out1.csv"), DataFrame)
            @test df.SeqName == ["seq9", "seq5", "seq1"]       # order of first appearance
            @test df.A == [4, 5, 16] && df.Sequence == ["AAA", "GGG", "TTT"]
        finally
            rm(td; recursive=true, force=true)
        end
    end

    @testset "merge_taxonomy_counts output ignores input order, including unnumbered names" begin
        db = SV.DatabaseMeta("silva", ["Domain", "Phylum"], "silva", Dict{String,Any}[], Set{String}())
        hits = ["seq10\tBacteria;Firmicutes\t95.0", "seq2\tBacteria;Proteobacteria\t88.0",
                "ASVb\tBacteria;Actino\t91.0", "ASVa\tEukaryota;Chloro\t99.0", "otu3\tBacteria;Bact\t90.0"]
        counts = ["seq10,5,1", "seq2,3,0", "ASVb,1,1", "ASVa,9,9", "otu3,0,4", "zzz_counts_only,2,2"]
        results = String[]
        for perm in (1:5, [5, 3, 1, 4, 2], [2, 4, 5, 1, 3])
            v = tempname() * ".tsv"; c = tempname() * ".csv"
            write(v, join(hits[collect(perm)], "\n") * "\n")
            cperm = vcat(collect(perm), 6)
            write(c, "SeqName,s1,s2\n" * join(counts[reverse(cperm)], "\n") * "\n")
            df = SV.TaxonomyTableTools.merge_taxonomy_counts(v, c, db)
            io = IOBuffer(); CSV.write(io, df); push!(results, String(take!(io)))
            rm(v); rm(c)
        end
        @test allequal(results)
        df = CSV.read(IOBuffer(results[1]), DataFrame)
        @test df.SeqName == ["seq2", "otu3", "seq10", "ASVa", "ASVb", "zzz_counts_only"]
    end

    @testset "composition summary is byte-stable and ordered by reads" begin
        root = mktempdir()
        try
            SV.ServerState.set_root!(root)
            mkpath(joinpath(root, "config"))
            write(joinpath(root, "config", "composition.yml"), """
            filters:
              euk: { filters: [ { column: Genus, type: include, values: [Blastocystis, Giardia, Entamoeba] } ] }
            sets:
              s: { label: S, categories: [ { name: Euk, filter: euk }, { name: Other } ] }
            """)
            bodies = String[]
            for shuffled in (false, true)
                run = shuffled ? "shuf" : "sorted"
                mkpath(joinpath(root, "data", "st", run))
                md = joinpath(root, "projects", "st", run, "merged"); mkpath(md)
                mem, mcon = _det_db(; threads=last(_DET_THREADS), preserve=false, shuffled)
                DBInterface.execute(mcon, "ATTACH '$(joinpath(md, "results.duckdb"))' AS f")
                DBInterface.execute(mcon, "CREATE TABLE f.merged AS SELECT * FROM merged")
                DBInterface.close!(mcon); close(mem)
                for _ in 1:2
                    push!(bodies, String(SV._composition_summary("st", run, "s", nothing).body))
                end
            end
            @test allequal(bodies)
            cats = JSON3.read(bodies[1]).categories
            reads = [c.reads for c in values(cats)]
            @test issorted(reads; rev=true)
        finally
            rm(root; recursive=true, force=true)
        end
    end

    @testset "R ordination and PERMANOVA are seeded, and immune to session RNG state" begin
        if !SV.r_available()
            @info "R/vegan unavailable - skipping NMDS/PERMANOVA determinism"
        else
            rng = MersenneTwister(11)
            mat = Float64.(rand(rng, 0:50, 12, 60))
            meta = DataFrame(sample=["s$i" for i in 1:12], grp=repeat(["a", "b", "c"], 4))

            base_nmds = SV.run_nmds(mat; seed=123)
            base_perm = SV.run_permanova(mat, meta; seed=123)
            @test !any(isnan, base_nmds[1])
            @test base_perm isa NamedTuple && haskey(base_perm, :p_value)
            @test [g.group for g in base_perm.dispersion.groups] == ["a", "b", "c"]
            @test base_perm.dispersion.df == [2, 9]
            @test 0 < base_perm.dispersion.p_value <= 1

            @test SV.run_nmds(mat; seed=123) == base_nmds
            @test SV.run_permanova(mat, meta; seed=123) == base_perm

            # Perturb the embedded session the way another caller or a user's
            # .Rprofile could: a different RNG kind, an advanced stream, and a
            # parallel default. Results must not move.
            SV.MetaManifold.RRuntime.with_r_lock() do
                SV.Analysis.RCall.reval("""RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rounding");
                                  invisible(runif(1000)); options(mc.cores = 4)""")
            end
            try
                @test SV.run_nmds(mat; seed=123) == base_nmds
                @test SV.run_permanova(mat, meta; seed=123) == base_perm
            finally
                SV.MetaManifold.RRuntime.with_r_lock() do
                    SV.Analysis.RCall.reval("""suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"));
                                      options(mc.cores = NULL)""")
                end
            end

            # The seed is load-bearing: the permutation p-value is a Monte Carlo
            # estimate, so some other seed must move it.
            @test any(s -> SV.run_permanova(mat, meta; seed=s).p_value != base_perm.p_value, 1:10)
        end
    end
end

# Testsets above point the server at temporary roots.
SV.ServerState._root[] = _saved_root
