# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Phylogenetic placement: FASTA intake, settings, the command lines each step
# builds, where each step runs, the QC summaries, and the routes.
using MetaManifold
using JSON3, YAML, DuckDB, DBInterface, Random

const _PH = MetaManifold.Phylogeny
const _PSV = MetaManifold.Server

_ph_request(method, path, body="", headers=["Content-Type" => "application/json"]) =
    _PSV.Oxygen.internalrequest(_PSV.HTTP.Request(method, path, headers, body))
_ph_json(r) = JSON3.read(String(r.body))

_ph_factory() = YAML.load_file(joinpath(@__DIR__, "..", "..", "config", "defaults", "pipeline.yml"))

# Sequences related by a simple tree: two clades of four, and query fragments from each.
function _ph_sequences(rng)
    rand_seq(n) = join(rand(rng, ['A', 'C', 'G', 'T'], n))
    mutate(s, k) = (c = collect(s); for i in rand(rng, 1:length(c), k); c[i] = rand(rng, ['A', 'C', 'G', 'T']); end; join(c))
    root = rand_seq(400)
    a, b = mutate(root, 60), mutate(root, 60)
    refs = [("cladeA_$i" => mutate(a, 12)) for i in 1:4]
    append!(refs, [("cladeB_$i" => mutate(b, 12)) for i in 1:4])
    queries = ["qA" => mutate(a, 6)[101:300], "qB" => mutate(b, 6)[51:250]]
    refs, queries
end

@testset "Phylogeny" begin

    @testset "FASTA names are made safe and sequences ungapped" begin
        recs = _PH.parse_fasta(">Tritrichomonas foetus (AF1)\nACGT-\nac.gt\n>b:c,d\nNNNN\n")
        @test recs == ["Tritrichomonas_foetus__AF1_" => "ACGTACGT", "b_c_d" => "NNNN"]
        @test_throws ErrorException _PH.parse_fasta(">a b\nACGT\n>a_b\nACGT\n")
        @test_throws ErrorException _PH.parse_fasta(">a\nACGTX\n")
        @test_throws ErrorException _PH.parse_fasta(">a\n\n>b\nACGT\n")
        @test_throws ErrorException _PH.parse_fasta("ACGT\n")
        @test_throws ErrorException _PH.parse_fasta("")
    end

    @testset "overrides lie over the configured settings" begin
        cfg = _ph_factory()
        s = _PH.placement_settings(cfg, Dict("reference" => Dict("trim" => Dict("gap_threshold" => 0.05))))
        @test s["reference"]["trim"]["gap_threshold"] == 0.05
        @test s["reference"]["trim"]["method"] == "manual"
        @test s["placement"]["trim"]["gap_threshold"] == 0.01
        bad(ov) = @test_throws ErrorException _PH.placement_settings(cfg, ov)
        bad(Dict("reference" => Dict("tree" => Dict("bootstrap" => "ultrafast", "replicates" => 100))))
        bad(Dict("placement" => Dict("accumulate" => Dict("threshold" => 0.3))))
        bad(Dict("reference" => Dict("trim" => Dict("method" => "tight"))))
        bad(Dict("reference" => Dict("trim" => Dict("residue_overlap" => 0.75))))
        bad(Dict("placement" => Dict("trim" => Dict("conservation" => 150))))
    end

    @testset "trimAl flags follow the trim method" begin
        @test _PH.trim_args(Dict("method" => "manual", "gap_threshold" => 0.3)) == "-gt 0.3"
        @test _PH.trim_args(Dict("method" => "manual", "gap_threshold" => 0.8, "conservation" => 60,
                                 "similarity_threshold" => 0.001)) == "-gt 0.8 -cons 60 -st 0.001"
        @test _PH.trim_args(Dict("method" => "gappyout", "gap_threshold" => 0.3)) == "-gappyout"
        @test _PH.trim_args(Dict("method" => "strict", "residue_overlap" => 0.75,
                                 "sequence_overlap" => 80)) == "-strict -resoverlap 0.75 -seqoverlap 80"
    end

    @testset "each step builds the commands of the workflow" begin
        s = _PH.placement_settings(_ph_factory())
        bins = Dict("mafft" => "mafft", "trimal" => "trimal", "iqtree" => "iqtree",
                    "raxml" => "raxmlHPC-PTHREADS-SSE3", "gappa" => "gappa")
        cmd(steps, i; s=s) = _PH.step_commands(steps[i], s, bins, 8, 123; path=identity, workdir="/stage/x/")
        R, P = _PH.REFERENCE_STEPS, _PH.PLACEMENT_STEPS
        @test cmd(R, 1) == ["mafft --thread 8 --maxiterate 1000 --localpair references.fasta > reference.aln.fasta"]
        @test cmd(R, 2) == ["trimal -in reference.aln.fasta -out reference.trim.fasta -fasta -gt 0.3 " *
                            "-colnumbering > reference.columns"]
        @test cmd(R, 3) == ["iqtree -s reference.trim.fasta -m MFP -b 100 -nt 8 -seed 123 -pre reference -redo"]
        @test cmd(P, 1) == ["mafft --thread 8 --auto --addfragments queries.fasta reference.trim.fasta > combined.aln.fasta"]
        @test cmd(P, 2) == ["trimal -in combined.aln.fasta -out combined.trim.fasta -fasta -gt 0.01 " *
                            "-colnumbering > combined.columns"]
        @test cmd(P, 3) == ["raxmlHPC-PTHREADS-SSE3 -f v -m GTRCATI -G 0.2 -n epa -s combined.trim.fasta " *
                            "-t reference.treefile -T 8 -w /stage/x/"]
        @test cmd(P, 4) == ["gappa edit accumulate --jplace-path placement.jplace --threshold 0.8 " *
                            "--out-dir accumulate --allow-file-overwriting --threads 8"]
        ov = _PH.placement_settings(_ph_factory(), Dict(
            "reference" => Dict("tree" => Dict("bootstrap" => "ultrafast", "replicates" => 1000)),
            "placement" => Dict("align" => Dict("strategy" => "localpair", "maxiterate" => 1000),
                                "place" => Dict("heuristic" => nothing))))
        @test occursin("-bb 1000", only(cmd(R, 3; s=ov)))
        @test startswith(only(cmd(P, 1; s=ov)), "mafft --thread 8 --maxiterate 1000 --localpair --addfragments")
        @test !occursin("-G", only(cmd(P, 3; s=ov)))
        # RAxML's PTHREADS builds refuse a single thread.
        @test occursin("-T 2", only(_PH.step_commands(P[3], s, bins, 1, 1; path=identity, workdir="/w/")))
        # Extra flags pass the same gate as every other stage's.
        @test_throws ErrorException _PH.step_commands(R[3],
            _PH.placement_settings(_ph_factory(), Dict("reference" => Dict("tree" => Dict("optional_args" => "-x; rm")))),
            bins, 1, 1; path=identity)
    end

    @testset "a step goes to the server only when remote.stages lists it" begin
        cfg = _ph_factory()
        @test isnothing(_PH.remote_step_target(cfg, "phylogeny_tree"; threads=4))
        cfg["remote"] = Dict("host" => "me@server", "staging_dir" => "/scratch/stage",
                             "stages" => ["phylogeny_tree", "phylogeny_place"],
                             "tools" => Dict("iqtree" => "/opt/iqtree/bin/iqtree2", "raxml" => nothing))
        t = _PH.remote_step_target(cfg, "phylogeny_tree"; threads=4)
        @test t.host == "me@server"
        @test t.threads == 4
        @test t.tools["iqtree"] == "/opt/iqtree/bin/iqtree2"
        @test _PH.remote_step_target(cfg, "phylogeny_place"; threads=4).tools["raxml"] == "raxmlHPC-PTHREADS-SSE3"
        @test isnothing(_PH.remote_step_target(cfg, "phylogeny_align"; threads=4))
        @test isnothing(_PH.remote_step_target(cfg, nothing; threads=4))
        cfg["remote"]["threads"] = 24
        @test _PH.remote_step_target(cfg, "phylogeny_tree"; threads=4).threads == 24
        cfg["remote"]["threads"] = true
        @test _PH.remote_step_target(cfg, "phylogeny_tree"; threads=4).threads == 4
        cfg["remote"]["tools"]["iqtree"] = "iqtree; rm -rf ~"
        @test_throws ErrorException _PH.remote_step_target(cfg, "phylogeny_tree"; threads=4)
        # Trimming and accumulation never leave this machine.
        @test all(isnothing(s.remote) for s in (_PH.REFERENCE_STEPS[2], _PH.PLACEMENT_STEPS[2], _PH.PLACEMENT_STEPS[4]))
    end

    @testset "validation knows the phylogeny keys" begin
        @test "phylogeny_add" in Validation.REMOTE_STAGES
        cfg = _ph_factory()
        errs(c) = (e = Validation.ValidationError[]; Validation._validate_pipeline_cfg(e, c, "t"); e)
        @test isempty(errs(cfg))
        bad = deepcopy(cfg)
        bad["remote"]["tools"]["mafft"] = "-oProxyCommand=x"
        bad["phylogeny"]["placement"]["trim"]["gap_threshold"] = 2
        bad["phylogeny"]["reference"]["tree"]["model"] = "GTR; ls"
        msgs = [e.message for e in errs(bad)]
        @test any(occursin("remote.tools.mafft", m) for m in msgs)
        @test any(occursin("phylogeny.placement.trim.gap_threshold", m) for m in msgs)
        @test any(occursin("phylogeny.reference.tree.model", m) for m in msgs)
    end

    @testset "QC reads alignments, trees and placements" begin
        mktempdir() do d
            aln = joinpath(d, "a.fasta")
            write(aln, ">r1\nACGT-A\n>r2\nAC-T-A\n>q1\n--GT--\n>q2\n----C-\n")
            trimmed = joinpath(d, "t.fasta")
            write(trimmed, ">r1\nACTA\n>r2\nACTA\n>q1\n--T-\n")
            qc = _PH.alignment_qc(aln; trimmed, columns=[0, 1, 3, 5], queries=Set(["q1", "q2"]))
            @test qc["sequences"] == 4
            @test qc["columns"] == 6
            @test qc["occupancy"] == [0.5, 0.5, 0.5, 0.75, 0.25, 0.5]
            @test qc["occupancy_queries"] == [0.0, 0.0, 0.5, 0.5, 0.5, 0.0]
            @test qc["removed"] == ["q2"]
            q1 = only(filter(p -> p["name"] == "q1", qc["per_sequence"]))
            @test (q1["residues"], q1["kept_residues"], q1["span"], q1["query"]) == (2, 1, [2, 3], true)
            write(joinpath(d, "c.txt"), "#ColumnsMap\t0, 1, 3, 5\n")
            @test _PH.read_columns(joinpath(d, "c.txt")) == [0, 1, 3, 5]
        end
        @test _PH._supports("((a:1,b:1)95:0.1,(c:1,d:1)70/88:0.2,e:1)100;") == [95.0, 70.0, 100.0]
        mktempdir() do d
            jp = joinpath(d, "p.jplace")
            write(jp, """{"version":3,"tree":"((a:1{0},b:1{1}):1{2},c:1{3});","fields":["edge_num","likelihood","like_weight_ratio","distal_length","pendant_length"],
                "placements":[{"p":[[0,-10,0.7,0.1,0.1],[1,-11,0.3,0.1,0.1]],"n":["q1"]},{"p":[[3,-9,1.0,0.1,0.1]],"nm":[["q2",1]]}]}""")
            q = _PH._jplace_queries(jp)
            @test q["q1"] == Dict("placements" => 2, "best_lwr" => 0.7, "edge" => 0)
            @test q["q2"]["edge"] == 3
        end
    end

    @testset "library and placement routes" begin
        tmp = mktempdir()
        old_root = _PSV.ServerState._root[]
        _PSV.ServerState.set_root!(tmp)
        mkpath(joinpath(tmp, "config", "defaults"))
        cp(joinpath(@__DIR__, "..", "..", "config", "defaults", "pipeline.yml"),
           joinpath(tmp, "config", "defaults", "pipeline.yml"))
        mkpath(joinpath(tmp, "data", "StudyP", "run1"))
        touch(joinpath(tmp, "data", "StudyP", "run1", "s1_R1.fastq.gz"))
        merged = joinpath(tmp, "projects", "StudyP", "run1", "merged")
        mkpath(merged)
        db = DBInterface.connect(DuckDB.DB, joinpath(merged, "results.duckdb"))
        DBInterface.execute(db, """CREATE TABLE merged AS SELECT * FROM (VALUES
            ('seq1', 'ACGTACGT', 'Parabasalia', 5, 0),
            ('seq2', 'ACGTTTGT', 'Parabasalia', 0, 3),
            ('seq3', 'GGGTACGT', 'Fornicata',   2, 1)) t(SeqName, Sequence, Class, Caecum_1, Colon_1)""")
        DBInterface.close!(db)
        try
            lib = "/api/v1/reference-trees"
            @test isempty(_ph_json(_ph_request("GET", lib)))
            @test _ph_request("POST", lib, JSON3.write(Dict("name" => ""))).status == 400
            ref = _ph_json(_ph_request("POST", lib, JSON3.write(Dict("name" => "Parabasalia refs")))).id
            r = _ph_request("PUT", "$lib/$ref/fasta", ">Ref one\nACGT\n>r2\nACGA\n", ["Content-Type" => "text/plain"])
            @test (_ph_json(r).count, _ph_json(r).renamed) == (2, 1)
            @test String(_ph_request("GET", "$lib/$ref/fasta").body) == ">Ref_one\nACGT\n>r2\nACGA\n"
            @test _ph_request("PUT", "$lib/$ref", JSON3.write(Dict("overrides" =>
                Dict("reference" => Dict("trim" => Dict("method" => "gappyout")))))).status == 200
            @test _ph_request("PUT", "$lib/$ref", JSON3.write(Dict("overrides" =>
                Dict("placement" => Dict("trim" => Dict("gap_threshold" => 0.1)))))).status == 400
            d = _ph_json(_ph_request("GET", "$lib/$ref"))
            @test d.settings.reference.trim.method == "gappyout"
            @test d.inherited.reference.trim.method == "manual"
            @test d.state == "new"
            @test _ph_request("GET", "$lib/$ref/qc/align").status == 404
            @test _ph_request("GET", "$lib/$ref/tree").status == 404
            mkpath(joinpath(tmp, "reference_trees", ref, "tree"))
            write(joinpath(tmp, "reference_trees", ref, "tree", "reference.treefile"), "((a:1,b:1)90:1,c:1,d:1);")
            t = _ph_json(_ph_request("GET", "$lib/$ref/tree"))
            @test (t.file, t.format, t.view) == ("Parabasalia_refs.treefile", "newick", nothing)
            @test _ph_request("PUT", "$lib/$ref/tree/view", JSON3.write(Dict("version" => 1))).status == 200
            @test _ph_json(_ph_request("GET", "$lib/$ref/tree")).view.version == 1
            @test _ph_request("PUT", "$lib/$ref/tree/view", "[1]").status == 400
            rm(joinpath(tmp, "reference_trees", ref, "tree"); recursive=true)
            @test _ph_request("POST", "$lib/$ref/trim-preview", JSON3.write(Dict("method" => "manual"))).status == 400

            base = "/api/v1/studies/StudyP/placements"
            @test _ph_request("POST", base, JSON3.write(Dict("name" => "P", "reference" => "00000000"))).status == 400
            id = _ph_json(_ph_request("POST", base, JSON3.write(Dict("name" => "Parabasalia", "reference" => ref)))).id
            @test only(_ph_json(_ph_request("GET", base))).state == "new"
            @test _ph_request("DELETE", "$lib/$ref").status == 409

            q = Dict("source" => "taxon", "runs" => [Dict("run" => "run1")], "table" => "merged",
                     "rank" => "Class", "values" => ["Parabasalia"], "min_reads" => 1)
            @test _ph_json(_ph_request("POST", "$base/preview", JSON3.write(q))).count == 2
            q["runs"] = [Dict("run" => "run1", "subgroups" => ["Caecum"])]
            p = _ph_json(_ph_request("POST", "$base/preview", JSON3.write(q)))
            @test p.count == 1
            @test only(p.per_run).subgroups == ["Caecum"]
            q["runs"] = [Dict("run" => "run1", "subgroups" => ["Ileum"])]
            @test _ph_request("POST", "$base/preview", JSON3.write(q)).status == 400
            q["runs"] = [Dict("run" => "run1", "subgroups" => ["Colon"])]
            @test first.(_PSV._write_taxon_queries!("StudyP", id, q).records) == ["run1_seq2"]

            @test _ph_request("PUT", "$base/$id", JSON3.write(Dict("queries" => q,
                "overrides" => Dict("placement" => Dict("trim" => Dict("gap_threshold" => 0.05)))))).status == 200
            @test _ph_request("PUT", "$base/$id", JSON3.write(Dict(
                "overrides" => Dict("reference" => Dict("tree" => Dict("replicates" => 10)))))).status == 400
            @test _ph_request("POST", "$base/$id/run").status == 400
            d = _ph_json(_ph_request("GET", "$base/$id"))
            @test d.settings.placement.trim.gap_threshold == 0.05
            @test all(isnothing, values(d.remote))

            @test _ph_request("GET", "$base/nothex12").status == 400
            @test _ph_request("DELETE", "$base/$id").status == 200
            @test _ph_request("DELETE", "$lib/$ref").status == 200
            @test isempty(_ph_json(_ph_request("GET", lib)))
        finally
            _PSV.ServerState.set_root!(old_root)
            rm(tmp; recursive=true, force=true)
        end
    end

    # Runs the real tools when this machine has them all.
    have_tools = all(t -> !isnothing(Sys.which(MetaManifold.Tools.tool_bin(t))),
                     ("mafft", "trimal", "iqtree", "raxml", "gappa"))
    @testset "a reference tree and a placement run end to end" begin
        if !have_tools
            @info "Phylogeny: skipping the end-to-end run; not every tool is installed"
        else
            mktempdir() do d
                refs, queries = _ph_sequences(MersenneTwister(7))
                rdir, pdir = joinpath(d, "ref"), joinpath(d, "pl")
                rfiles = _PH.reference_files(rdir)
                pfiles = _PH.placement_files(pdir, rdir)
                _PH.write_fasta(rfiles["references.fasta"], refs)
                _PH.write_fasta(pfiles["queries.fasta"], queries)
                cfg = _ph_factory()
                ov = Dict("reference" => Dict("tree" => Dict("replicates" => 10, "model" => "GTR+G")))
                _PH.check_reference(rfiles)
                _PH.run_workflow(_PH.REFERENCE_STEPS, rdir, rfiles, cfg; overrides=ov)
                @test _PH.read_status(rdir)["state"] == "done"
                trimqc = JSON3.read(read(joinpath(rdir, "qc", "trim.json"), String))
                @test trimqc.sequences == 8
                @test !isempty(trimqc.kept_columns)
                treeqc = JSON3.read(read(joinpath(rdir, "qc", "tree.json"), String))
                @test !isnothing(treeqc.model)
                @test length(treeqc.supports) >= 4

                preview = _PH.trim_preview(_PH.REFERENCE_STEPS, rdir, rfiles, Dict("method" => "gappyout"))
                @test preview["sequences"] == 8
                @test _PH.read_status(rdir)["state"] == "done"

                _PH.check_placement(pfiles)
                _PH.run_workflow(_PH.PLACEMENT_STEPS, pdir, pfiles, cfg)
                plqc = JSON3.read(read(joinpath(pdir, "qc", "place.json"), String))
                @test plqc.placed == 2
                @test isfile(pfiles["accumulated.jplace"])

                # A new trim reruns trim and what follows it, and leaves the MAFFT alignment alone.
                aln_time = mtime(pfiles["combined.aln.fasta"])
                _PH.run_workflow(_PH.PLACEMENT_STEPS, pdir, pfiles, cfg;
                                 overrides=Dict("placement" => Dict("trim" => Dict("method" => "noallgaps"))))
                st = _PH.read_status(pdir)["steps"]
                @test st["align"]["state"] == "current"
                @test st["trim"]["state"] == "done"
                @test mtime(pfiles["combined.aln.fasta"]) == aln_time
            end
        end
    end
end
