# Sample total-read filtering.
#
# A sample is a column in the merged table, so a floor on a sample's total reads
# drops columns, not rows. These tests pin the helper semantics and then drive
# every HTTP surface that honours the filter through Oxygen's in-process router,
# on a fixture built around the canonical use:
#
#   exclude the Contaminant category  AND  Genus_boot >= 80  AND  sample reads >= 50
#
# Fixture (rows inserted out of SeqName order on purpose):
#
#   SeqName  Genus         Genus_boot  contamination  Srich  Slow  Scontam  Sdeep
#   seq1     Blastocystis  95          Retained          30    40        5    100
#   seq2     Giardia       60          Retained          20     0        0     10
#   seq3     Escherichia   99          Contaminant      500     1      900      0
#   seq4     Blastocystis  85          Retained          25     5       10      0
#
# Rows surviving (not Contaminant, boot >= 80): seq1, seq4.
#   per-sample reads over those rows ("filtered" basis): Srich 55, Slow 45, Scontam 15, Sdeep 100
#   per-sample library size ("raw" basis):               Srich 575, Slow 46, Scontam 915, Sdeep 110

using MetaManifold
SV = MetaManifold.Server

const _SR_SAMPLES = ["Srich", "Slow", "Scontam", "Sdeep"]

const _SR_LIBRARY = """
filters:
  whitelist:
    filters:
      - { column: Genus, type: include, values: [Blastocystis, Giardia] }
sets:
  contamination:
    label: Contamination
    categories:
      - { name: Retained, filter: whitelist }
      - { name: Contaminant }
"""

# Build <root>/projects/<study>/<run>/merged/results.duckdb plus the data-side
# run directory the route guards look for.
function _sr_make_run(root::String, study::String, run::String)
    data_run = joinpath(root, "data", study, run)
    mkpath(data_run)
    touch(joinpath(data_run, "Srich_R1.fastq.gz")); touch(joinpath(data_run, "Srich_R2.fastq.gz"))
    # No inherited figure exclusions: the request's filters are the only ones.
    write(joinpath(data_run, "pipeline.yml"), "analysis:\n  exclude_categories: []\n")

    merge_dir = joinpath(root, "projects", study, run, "merged")
    mkpath(merge_dir)
    db = DuckDB.DB(joinpath(merge_dir, "results.duckdb"))
    con = DBInterface.connect(db)
    try
        DBInterface.execute(con, """
            CREATE TABLE merged (
                "SeqName" VARCHAR, "sequence" VARCHAR, "Pident" DOUBLE,
                "Domain" VARCHAR, "Genus" VARCHAR, "Genus_boot" DOUBLE,
                "Srich" BIGINT, "Slow" BIGINT, "Scontam" BIGINT, "Sdeep" BIGINT)
        """)
        DBInterface.execute(con, """
            INSERT INTO merged VALUES
              ('seq4','ACGTTT',99.0,'Eukaryota','Blastocystis',85.0, 25, 5, 10,   0),
              ('seq2','ACGTCC',97.0,'Eukaryota','Giardia',     60.0, 20, 0,  0,  10),
              ('seq1','ACGTAA',98.0,'Eukaryota','Blastocystis',95.0, 30,40,  5, 100),
              ('seq3','ACGTGG',99.5,'Bacteria', 'Escherichia', 99.0,500, 1,900,   0)
        """)
        lib = SV.CompositionLibrary.load(joinpath(root, "config", "composition.yml"))
        SV.Categories.write_category_columns!(con, "merged", "VSEARCH", ["contamination"];
                                              library=lib)
    finally
        DBInterface.close!(con); close(db)
    end
    merge_dir
end

function _sr_fixture(f::Function)
    root = mktempdir()
    old_root = SV.ServerState._root[]
    try
        SV.ServerState.set_root!(root)
        mkpath(joinpath(root, "config"))
        write(joinpath(root, "config", "composition.yml"), _SR_LIBRARY)
        _sr_make_run(root, "studyS", "runA")
        _sr_make_run(root, "studyS", "runB")
        f(root)
    finally
        SV.ServerState._root[] = old_root
        rm(root; recursive=true, force=true)
    end
end

_sr_body(d) = JSON3.read(JSON3.write(d))

function _sr_post(path::String, body)
    req = SV.HTTP.Request("POST", path, ["Content-Type" => "application/json"], JSON3.write(body))
    SV.Oxygen.internalrequest(req)
end

# The user's example, expressed the way the Tables UI sends it.
_sr_filters(sample_reads) = Dict(
    "Category__contamination" => Dict("exclude" => ["Contaminant"]),
    "Genus_boot"              => Dict("min" => 80),
    SV.SAMPLE_READS_FILTER_KEY => sample_reads,
)

@testset "Sample total-read filter" begin

    @testset "bounds parsing" begin
        b = SV._sample_read_bounds(Dict("sample_min_reads" => "50"))
        @test b.min == 50.0 && isnothing(b.max) && b.basis == "filtered"
        @test SV._has_sample_read_bounds(b)

        b = SV._sample_read_bounds(Dict("sample_max_reads" => "10.5", "sample_reads_basis" => "raw"))
        @test isnothing(b.min) && b.max == 10.5 && b.basis == "raw"

        # An unknown basis falls back to the default rather than erroring, and a
        # non-numeric bound is ignored rather than filtering everything out.
        b = SV._sample_read_bounds(Dict("sample_min_reads" => "abc", "sample_reads_basis" => "bogus"))
        @test isnothing(b.min) && b.basis == "filtered"
        @test !SV._has_sample_read_bounds(b)
        @test !SV._has_sample_read_bounds(SV._sample_read_bounds(Dict{String,String}()))
    end

    @testset "request body: reserved colFilters key and explicit sampleReads" begin
        p = SV._body_filter_params(_sr_body(Dict("colFilters" => Dict(
            SV.SAMPLE_READS_FILTER_KEY => Dict("min" => 50, "basis" => "raw"),
            "Genus_boot" => Dict("min" => 80)))))
        @test p["sample_min_reads"] == "50"
        @test p["sample_reads_basis"] == "raw"
        @test p["col_min.Genus_boot"] == "80"
        # The reserved key must never become a row filter on a column of that name.
        @test !any(startswith(k, "col") && occursin(SV.SAMPLE_READS_FILTER_KEY, k) for k in keys(p))

        # Explicit top-level sampleReads wins over the reserved entry.
        p = SV._body_filter_params(_sr_body(Dict(
            "sampleReads" => Dict("min" => 7),
            "colFilters"  => Dict(SV.SAMPLE_READS_FILTER_KEY => Dict("min" => 50, "max" => 99)))))
        @test p["sample_min_reads"] == "7"
        @test p["sample_max_reads"] == "99"   # not set explicitly, so taken from the entry
    end

    @testset "_retain_sample_columns semantics" begin
        _sr_fixture() do root
            merge_dir = joinpath(root, "projects", "studyS", "runA", "merged")
            SV.with_results_db(merge_dir) do con
                cols = SV._duckdb_columns(con, "merged")
                scols = SV._sample_count_columns(con, "merged")
                @test scols == _SR_SAMPLES     # ordinal order, bootstraps and Pident excluded

                params = SV._body_filter_params(_sr_body(Dict("colFilters" => _sr_filters(Dict("min" => 50)))))
                where, wp = SV._build_where(params, cols)
                totals = SV._sample_read_totals(con, "merged", scols, where, wp)
                @test totals == Dict("Srich" => 55.0, "Slow" => 45.0, "Scontam" => 15.0, "Sdeep" => 100.0)
                @test SV._sample_read_totals(con, "merged", scols) ==
                      Dict("Srich" => 575.0, "Slow" => 46.0, "Scontam" => 915.0, "Sdeep" => 110.0)

                keep(sr) = SV._retain_sample_columns(con, "merged", scols,
                    SV._body_filter_params(_sr_body(Dict("colFilters" => _sr_filters(sr)))), where, wp)
                @test keep(Dict("min" => 50))                     == ["Srich", "Sdeep"]
                @test keep(Dict("min" => 50, "basis" => "raw"))   == ["Srich", "Scontam", "Sdeep"]
                @test keep(Dict("max" => 60))                     == ["Srich", "Slow", "Scontam"]
                @test keep(Dict("min" => 15, "max" => 55))        == ["Srich", "Slow", "Scontam"]  # inclusive
                @test keep(Dict("min" => 1000))                   == String[]
                # No bound: the input comes back untouched, with no query issued.
                @test SV._retain_sample_columns(con, "merged", scols, Dict{String,String}(), where, wp) === scols
            end
        end
    end

    @testset "Tables query route" begin
        _sr_fixture() do root
            url = "/api/v1/studies/studyS/runs/runA/results/tables/merged/query"
            r = _sr_post(url, Dict("page" => 1, "perPage" => 100,
                                   "colFilters" => _sr_filters(Dict("min" => 50))))
            @test r.status == 200
            d = JSON3.read(String(r.body))
            @test d.total == 2
            @test collect(d.sample_count_columns) == ["Srich", "Sdeep"]
            @test collect(d.excluded_samples) == ["Slow", "Scontam"]
            @test d.total_reads == 155              # 55 + 100
            @test d.total_reads_unfiltered == 1646  # every read in every sample
            @test !("Slow" in d.columns) && !("Scontam" in d.columns)
            @test "Genus_boot" in d.columns && "Category__contamination" in d.columns
            @test [row.SeqName for row in d.rows] == ["seq1", "seq4"]
            @test all(!haskey(row, :Slow) && !haskey(row, :Scontam) for row in d.rows)

            # Raw basis measures the library, so the contaminant-heavy sample survives.
            r = _sr_post(url, Dict("page" => 1, "perPage" => 100,
                                   "colFilters" => _sr_filters(Dict("min" => 50, "basis" => "raw"))))
            d = JSON3.read(String(r.body))
            @test collect(d.sample_count_columns) == ["Srich", "Scontam", "Sdeep"]
            @test d.total_reads == 170              # 55 + 15 + 100

            # No sample bound: nothing excluded, behaviour unchanged.
            r = _sr_post(url, Dict("page" => 1, "perPage" => 100))
            d = JSON3.read(String(r.body))
            @test collect(d.sample_count_columns) == _SR_SAMPLES
            @test isempty(d.excluded_samples)
            @test d.total == 4
        end
    end

    @testset "Composition query route (category column materialised on demand)" begin
        _sr_fixture() do root
            r = _sr_post("/api/v1/studies/studyS/runs/runA/composition/VSEARCH/query",
                         Dict("page" => 1, "perPage" => 100, "category_set" => "contamination",
                              "colFilters" => _sr_filters(Dict("min" => 50))))
            @test r.status == 200
            d = JSON3.read(String(r.body))
            @test collect(d.sample_count_columns) == ["Srich", "Sdeep"]
            @test d.total_reads == 155
        end
    end

    @testset "annotation-style select expression keeps its virtual column" begin
        _sr_fixture() do root
            merge_dir = joinpath(root, "projects", "studyS", "runA", "merged")
            SV.with_results_db(merge_dir) do con
                cols = SV._duckdb_columns(con, "merged")
                body = _sr_body(Dict("page" => 1, "perPage" => 10,
                                     "colFilters" => _sr_filters(Dict("min" => 50))))
                r = SV._duckdb_paginated_query(con, "merged", body;
                        select_expr = "*, '' AS \"Assignment\"",
                        response_columns = [cols; "Assignment"])
                d = JSON3.read(String(r.body))
                @test last(d.columns) == "Assignment"
                @test all(haskey(row, :Assignment) for row in d.rows)
                @test !("Slow" in d.columns)
            end
        end
    end

    @testset "Save and export routes" begin
        _sr_fixture() do root
            body = Dict("name" => "filtered_view", "colFilters" => _sr_filters(Dict("min" => 50)))
            r = _sr_post("/api/v1/studies/studyS/runs/runA/results/tables/merged/save", body)
            @test r.status == 200
            csv = CSV.read(joinpath(root, "projects", "studyS", "runA", "merged", "filtered_view.csv"), DataFrame)
            @test !("Slow" in names(csv)) && !("Scontam" in names(csv))
            @test "Srich" in names(csv) && "Sdeep" in names(csv)
            @test csv.SeqName == ["seq1", "seq4"]

            r = _sr_post("/api/v1/studies/studyS/runs/runA/results/tables/merged/export",
                         Dict("colFilters" => _sr_filters(Dict("min" => 50))))
            @test r.status == 200
            tmp = tempname() * ".xlsx"
            write(tmp, r.body)
            xf = SV.XLSX.readxlsx(tmp)
            sheet = xf[SV.XLSX.sheetnames(xf)[1]]
            header = [string(sheet[1, j]) for j in 1:size(sheet[:], 2)]
            @test !("Slow" in header) && "Sdeep" in header
            rm(tmp; force=true)

            r = _sr_post("/api/v1/studies/studyS/runs/runA/results/tables/nope/export", Dict())
            @test r.status == 404
        end
    end

    @testset "filter presets round-trip exclusions and the sample bound" begin
        _sr_fixture() do root
            filters = _sr_filters(Dict("min" => 50, "basis" => "raw"))
            r = _sr_post("/api/v1/filter-presets/sr_example",
                         Dict("filters" => filters, "description" => "example"))
            @test r.status == 200
            doc = YAML.load_file(joinpath(SV._presets_dir(), "sr_example.yml"))
            entries = doc["filters"]
            @test any(e -> get(e, "type", "") == "exclude" &&
                           e["column"] == "Category__contamination" &&
                           e["values"] == ["Contaminant"], entries)
            sr = only(filter(e -> get(e, "type", "") == "sample_reads", entries))
            @test sr["min"] == 50 && sr["basis"] == "raw"
            # Never saved as a bogus column filter.
            @test !any(e -> get(e, "column", "") == SV.SAMPLE_READS_FILTER_KEY, entries)

            r = _sr_post("/api/v1/studies/studyS/runs/runA/results/tables/merged/apply-preset",
                         Dict("preset" => "sr_example.yml"))
            @test r.status == 200
            d = JSON3.read(String(r.body))
            @test collect(d.filters.Category__contamination.exclude) == ["Contaminant"]
            @test d.filters.Genus_boot.min == 80
            restored = d.filters[Symbol(SV.SAMPLE_READS_FILTER_KEY)]
            @test restored.min == 50 && restored.basis == "raw"
            @test d.rows_before == 4
            @test d.rows_after == 2          # exclusion and bootstrap floor both applied

            # Re-sending the restored filters reproduces the original view.
            q = _sr_post("/api/v1/studies/studyS/runs/runA/results/tables/merged/query",
                         Dict("page" => 1, "perPage" => 100, "colFilters" => d.filters))
            @test collect(JSON3.read(String(q.body)).sample_count_columns) == ["Srich", "Scontam", "Sdeep"]
            rm(joinpath(SV._presets_dir(), "sr_example.yml"); force=true)
        end
    end

    @testset "Per-run chart and alpha routes" begin
        _sr_fixture() do root
            r = _sr_post("/api/v1/studies/studyS/runs/runA/analysis/chart",
                         Dict("table" => "merged", "tag" => "rank", "value" => "Genus",
                              "relative" => false, "colFilters" => _sr_filters(Dict("min" => 50))))
            @test r.status == 200
            fig = JSON3.read(String(r.body))
            xs = unique(vcat([collect(t.x) for t in fig.data]...))
            @test xs == ["Sdeep", "Srich"]

            r = _sr_post("/api/v1/studies/studyS/runs/runA/analysis/alpha",
                         Dict("table" => "merged", "colFilters" => _sr_filters(Dict("min" => 50))))
            @test r.status == 200
            fig = JSON3.read(String(r.body))
            @test collect(fig.data[1].x) == ["Sdeep", "Srich"]

            # A floor no sample reaches is a clear 400, not an empty figure.
            r = _sr_post("/api/v1/studies/studyS/runs/runA/analysis/chart",
                         Dict("table" => "merged", "tag" => "rank", "value" => "Genus",
                              "colFilters" => _sr_filters(Dict("min" => 10_000))))
            @test r.status == 400
        end
    end

    @testset "Composition summary route" begin
        _sr_fixture() do root
            r = _sr_post("/api/v1/studies/studyS/runs/runA/composition/summary",
                         Dict("category_set" => "contamination",
                              "colFilters" => Dict(SV.SAMPLE_READS_FILTER_KEY => Dict("min" => 50))))
            @test r.status == 200
            d = JSON3.read(String(r.body))
            # No row filters on this surface, so the bound reads library sizes.
            @test collect(d.samples) == ["Srich", "Scontam", "Sdeep"]
            @test d.total_reads == 575 + 915 + 110
            # Categories ordered by reads, descending.
            @test collect(keys(d.categories)) == [:Contaminant, :Retained]
        end
    end

    @testset "Cross-run chart, alpha, venn and ordination routes" begin
        _sr_fixture() do root
            runs = [Dict("run" => "runA"), Dict("run" => "runB")]
            filt = _sr_filters(Dict("min" => 50))

            r = _sr_post("/api/v1/studies/studyS/analysis/chart",
                         Dict("runs" => runs, "table" => "merged", "tag" => "rank",
                              "value" => "Genus", "relative" => false, "colFilters" => filt))
            @test r.status == 200
            fig = JSON3.read(String(r.body))
            # Cross-run bars pool each run into one column. Only the Blastocystis
            # rows survive the row filters, and only Srich + Sdeep survive the
            # sample floor: 55 + 100 = 155 reads per run. Without the floor, Slow
            # and Scontam would add 45 + 15 per run.
            @test [t.name for t in fig.data] == ["Blastocystis"]
            @test collect(fig.data[1].x) == ["runA", "runB"]
            @test collect(fig.data[1].y) == [155.0, 155.0]

            r = _sr_post("/api/v1/studies/studyS/analysis/alpha",
                         Dict("runs" => runs, "table" => "merged", "colFilters" => filt))
            @test r.status == 200
            s = String(r.body)
            @test !occursin("Slow", s) && occursin("Sdeep", s)

            # Venn works on presence: with only Srich and Sdeep kept, and rows
            # restricted to seq1/seq4, the single present genus is Blastocystis.
            r = _sr_post("/api/v1/studies/studyS/analysis/venn",
                         Dict("runs" => runs, "table" => "merged", "rank" => "Genus", "colFilters" => filt))
            @test r.status == 200
            @test !occursin("Giardia", String(r.body))

            if SV.r_available()
                r = _sr_post("/api/v1/studies/studyS/analysis/nmds",
                             Dict("runs" => runs, "table" => "merged", "colFilters" =>
                                  Dict(SV.SAMPLE_READS_FILTER_KEY => Dict("min" => 50))))
                # 3 samples per run survive (raw == filtered here: no row filters),
                # giving 6 points; the dropped sample must not be plotted.
                @test r.status == 200
                @test !occursin("Slow", String(r.body))
            else
                @info "R/vegan unavailable - skipping NMDS sample-filter route test"
            end
        end
    end
end
