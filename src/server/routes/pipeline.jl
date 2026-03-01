# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: pipeline execution
#
# Every job begins by calling write_run_config() for the target run(s) so that
# any config changes made in the UI are picked up before execution. The
# existing checkpoint/hash logic then skips stages whose inputs are unchanged.

using JSON3

const _r_lock    = ReentrantLock()
const _plot_lock = ReentrantLock()

## Stage execution helpers

function _project_ctx(study::String, run::String)
    projects = new_project(run;
                           data_dir     = joinpath(ServerState.data_dir(), study),
                           projects_dir = joinpath(ServerState.projects_dir(), study),
                           config_dir   = joinpath(dirname(ServerState.data_dir()), "config"))
    first(projects)
end

function _run_full_pipeline(study::String, run::String, db_config::String,
                             dbs::Dict, db_name::String)
    project  = _project_ctx(study, run)
    db_meta  = make_db_meta(db_config, db_name)
    trimmed  = cutadapt(project)

    dada2_task = Threads.@spawn begin
        asvs = lock(_r_lock) do
            dada2(project, trimmed, taxonomy_db=dbs["$(db_name)_dada2"])
        end
        cdhit(project, asvs)
    end
    swarm_task = Threads.@spawn swarm(project, trimmed)

    asvs = fetch(dada2_task)
    otus = fetch(swarm_task)

    asv_tax_task = Threads.@spawn vsearch(project, asvs, dbs["$(db_name)_vsearch"])
    otu_tax_task = isnothing(otus) ? nothing :
                  Threads.@spawn vsearch(project, otus, dbs["$(db_name)_vsearch"])

    asv_tax   = fetch(asv_tax_task)
    asv_merged = merge_taxa(project, asvs, asv_tax, db_meta)

    merged = if !isnothing(otus)
        otu_tax    = fetch(otu_tax_task)
        otu_merged = merge_taxa_otu(project, otus, otu_tax, db_meta)
        MergedTables(merge(asv_merged.tables, otu_merged.tables),
                     asv_merged.filter_order, asv_merged.filter_colours)
    else
        asv_merged
    end

    analyse_run(project, merged, asvs, db_meta; plot_lock=_plot_lock)
    merged
end

function _run_stage(study::String, run::String, stage::String,
                    db_config::String, dbs::Dict, db_name::String)
    project = _project_ctx(study, run)
    db_meta = make_db_meta(db_config, db_name)

    if stage == "cutadapt"
        cutadapt(project)
    elseif stage == "dada2_denoise"
        trimmed = TrimmedReads(joinpath(ServerState.projects_dir(), study, run, "cutadapt"))
        lock(_r_lock) do
            dada2_denoise(project, trimmed)
        end
    elseif stage == "dada2_classify"
        trimmed  = TrimmedReads(joinpath(ServerState.projects_dir(), study, run, "cutadapt"))
        denoised = DenoisedASVs(
            joinpath(ServerState.projects_dir(), study, run, "dada2", "Checkpoints", "ckpt_chimera.RData"),
            write_run_config(project),
            joinpath(ServerState.projects_dir(), study, run, "dada2"),
            trimmed.dir)
        lock(_r_lock) do
            dada2_classify(project, denoised; taxonomy_db=dbs["$(db_name)_dada2"])
        end
    elseif stage == "swarm"
        trimmed = TrimmedReads(joinpath(ServerState.projects_dir(), study, run, "cutadapt"))
        swarm(project, trimmed)
    elseif stage == "cdhit"
        asvs = _asvresult_from_disk(study, run)
        cdhit(project, asvs)
    elseif stage == "vsearch_asv"
        asvs = _asvresult_from_disk(study, run)
        vsearch(project, asvs, dbs["$(db_name)_vsearch"])
    elseif stage == "vsearch_otu"
        otus = _oturesult_from_disk(study, run)
        isnothing(otus) && error("swarm output not found - run swarm first")
        vsearch(project, otus, dbs["$(db_name)_vsearch"])
    elseif stage == "merge_taxa"
        asvs    = _asvresult_from_disk(study, run)
        asv_tax = TaxonomyHits(joinpath(ServerState.projects_dir(), study, run, "vsearch", "asv_hits.tsv"))
        merge_taxa(project, asvs, asv_tax, db_meta)
    elseif stage == "analyse_run"
        asvs   = _asvresult_from_disk(study, run)
        merged = _mergedtables_from_disk(study, run)
        analyse_run(project, merged, asvs, db_meta; plot_lock=_plot_lock)
    else
        error("Unknown stage: $stage")
    end
end

## Disk reconstruction helpers

function _asvresult_from_disk(study::String, run::String)
    t = joinpath(ServerState.projects_dir(), study, run, "dada2", "Tables")
    ASVResult(joinpath(t, "asvs.fasta"),
              joinpath(t, "seqtab_nochim.csv"),
              joinpath(t, "taxonomy.csv"))
end

function _oturesult_from_disk(study::String, run::String)
    p = joinpath(ServerState.projects_dir(), study, run, "swarm", "otu_table.csv")
    isfile(p) || return nothing
    f = joinpath(ServerState.projects_dir(), study, run, "swarm", "otus.fasta")
    OTUResult(f, p)
end

function _mergedtables_from_disk(study::String, run::String)
    merge_dir = joinpath(ServerState.projects_dir(), study, run, "merged")
    tables    = Dict{String,String}()
    isdir(merge_dir) || error("merge_taxa output not found - run merge_taxa first")
    for f in readdir(merge_dir; join=true)
        endswith(f, ".csv") || continue
        tables[splitext(basename(f))[1]] = f
    end
    MergedTables(tables, String[], Dict{String,String}())
end

## Routes

function _db_config_path()
    cfg = joinpath(dirname(ServerState.data_dir()), "config", "databases.yml")
    isfile(cfg) ? cfg : joinpath(dirname(ServerState.data_dir()), "config", "ci", "databases.yml")
end

function _load_dbs()
    path = _db_config_path()
    cfg  = YAML.load_file(path)
    db_name = get(get(cfg, "default", Dict()), "database", "pr2")
    dbs     = ensure_databases(path)
    db_name, dbs
end

@post "/api/v1/studies/{study}/pipeline" function(req, study::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    runs    = _run_names(study)
    db_cfg  = _db_config_path()

    job = submit_job!("pipeline"; study) do
        db_name, dbs = _load_dbs()
        merged_results = Vector{Any}(undef, length(runs))
        asvs_results   = Vector{Any}(undef, length(runs))

        Threads.@threads for (i, run) in collect(enumerate(runs))
            merged_results[i] = _run_full_pipeline(study, run, db_cfg, dbs, db_name)
        end

        valid = count(!isnothing, merged_results)
        if valid >= 2
            projects = [_project_ctx(study, r) for r in runs]
            db_meta  = make_db_meta(db_cfg, first(split(db_name, "_")))
            analyse_study(projects, merged_results, fill(db_meta, length(projects));
                          plot_lock=_plot_lock)
        end
    end

    json(job)
end

@post "/api/v1/studies/{study}/runs/{run}/pipeline" function(req, study::String, run::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    db_cfg = _db_config_path()

    job = submit_job!("pipeline"; study, run) do
        db_name, dbs = _load_dbs()
        _run_full_pipeline(study, run, db_cfg, dbs, db_name)
    end

    json(job)
end

@post "/api/v1/studies/{study}/runs/{run}/stages/{stage}" function(req, study::String,
                                                                        run::String,
                                                                        stage::String)
    study in _study_names() || return json_error(404, "study_not_found",
                                                     "Study '$study' not found")
    run in _run_names(study) || return json_error(404, "run_not_found",
                                                      "Run '$run' not found")
    stage in STAGES || return json_error(400, "unknown_stage", "Unknown stage: $stage")
    db_cfg = _db_config_path()

    job = submit_job!("stage"; study, run, stage) do
        db_name, dbs = _load_dbs()
        _run_stage(study, run, stage, db_cfg, dbs, db_name)
    end

    json(job)
end
