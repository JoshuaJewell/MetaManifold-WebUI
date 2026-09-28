# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Precompile workload: the requests a first visit to a study and a run makes,
# against a small temporary project, so the package cache holds the handlers'
# native code. Set `precompile_workload = false` for MetaManifold in
# LocalPreferences.toml to skip it while developing; the server's start-up
# warm-up then compiles the same routes at run time.
@setup_workload begin
    _defaults = joinpath(@__DIR__, "..", "..", "config", "defaults")
    _paths = let s = "/api/v1/studies/study", r = "/api/v1/studies/study/runs/run"
        ["/api/v1/studies", "/api/v1/jobs", s, "$s/config", "$s/config/overrides", "$s/runs",
         r, "$r/config", "$r/results/qc", "$r/results/dada2", "$r/results/tables",
         "$r/analysis/read-funnel"]
    end
    @compile_workload begin
        root = mktempdir()
        saved = ServerState._root[]
        try
            mkpath(joinpath(root, "config"))
            cp(_defaults, joinpath(root, "config", "defaults"))
            run_dir = joinpath(root, "data", "study", "run")
            mkpath(run_dir)
            foreach(r -> touch(joinpath(run_dir, "s1_L001_$(r)_001.fastq.gz")), ("R1", "R2"))
            merged = joinpath(root, "projects", "study", "run", "merged")
            mkpath(merged)
            db = DuckDB.DB(joinpath(merged, "results.duckdb"))
            con = DBInterface.connect(db)
            DBInterface.execute(con, "CREATE TABLE merged (SeqName VARCHAR, Sequence VARCHAR, Genus_dada2 VARCHAR, s1 BIGINT)")
            DBInterface.execute(con, "INSERT INTO merged VALUES ('seq1', 'ACGT', 'Giardia', 10)")
            DBInterface.close!(con)
            close(db)

            ServerState.set_root!(root)
            register_routes!()
            with_logger(NullLogger()) do
                for path in _paths
                    Oxygen.internalrequest(HTTP.Request("GET", path))
                end
                Oxygen.internalrequest(HTTP.Request("POST", "/api/v1/studies/study/runs/run/results/tables/merged/query",
                    ["Content-Type" => "application/json"], """{"page": 1, "page_size": 50}"""))
            end
        finally
            Oxygen.resetstate()
            ServerState._root[] = saved
            rm(root; recursive=true, force=true)
        end
    end
end
