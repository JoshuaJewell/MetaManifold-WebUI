#!/usr/bin/env julia
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

using Dates, SHA, JSON3, YAML
using MetaManifold.Provenance
include("biological_validation.jl")

function _tool_pin_matches(key::String, reported::String, pin::AbstractDict)
    expected = pin["version"]
    reported == expected && return true
    # The pinned VSEARCH binaries append their platform to the release number.
    # Retain that full identifier in provenance and accept only configured tags.
    key == "vsearch" || return false
    any(platform -> reported == expected * "_" * replace(String(platform), "-" => "_"),
        keys(get(pin, "archives", Dict())))
end

function _require_fresh_matrix(path::String, started::Real)
    isfinite(started) && started > 0 || error("invalid pipeline start time")
    isfile(path) || error("missing pipeline output $path; pipeline skipping is not validation")
    mtime(path) >= started || error("stale pipeline output $path; rerun the pipeline after recording its start time")
    filesize(path) > 0 || error("empty pipeline output $path")
    nothing
end

function collect_pipeline_evidence(output::String;
        started::Real=parse(Int, get(ENV, "STATISTICS_PIPELINE_STARTED_UNIX", "0")))
    isfinite(started) && 0 < started <= time() || error(
        "set STATISTICS_PIPELINE_STARTED_UNIX to the Unix timestamp recorded before starting the pipeline")
    root = normpath(joinpath(@__DIR__, "..", ".."))
    pins = YAML.load_file(joinpath(root, "config", "defaults", "tool_versions.yml"))
    string(VERSION) == pins["runtimes"]["julia"]["version"] || error("Julia does not match the repository pin")
    rversion = readchomp(`Rscript --vanilla -e 'cat(as.character(getRversion()))'`)
    rversion == pins["runtimes"]["r"]["version"] || error("R does not match the repository pin")
    dada2_version = readchomp(`Rscript --vanilla -e 'cat(as.character(packageVersion("dada2")))'`)
    lockfile = JSON3.read(read(joinpath(root, "renv.lock"), String))
    dada2_version == lockfile.Packages.dada2.Version || error("DADA2 does not match renv.lock")
    package_names = sort(String.(collect(keys(lockfile.Packages))))
    r_record = probe_r(; packages=package_names, renv_lock=joinpath(root, "renv.lock"))
    embedded_version = match(r"^R version ([0-9]+\.[0-9]+\.[0-9]+)(?: |$)", r_record.version)
    !isnothing(embedded_version) && embedded_version.captures[1] == rversion ||
        error("embedded R does not match the pinned Rscript runtime")
    for package in package_names
        r_record.packages[package] == lockfile.Packages[Symbol(package)].Version ||
            error("R package $package does not match renv.lock")
    end
    tool_records = Dict{String,Any}()
    for key in ("cutadapt", "fastqc", "multiqc", "vsearch", "swarm", "cd_hit_est")
        probe = TOOL_PROBES[key]
        record = probe_tool(probe; bin_resolver=_ -> probe.name)
        _tool_pin_matches(key, record.version, pins["tools"][key]) ||
            error("$key reports $(record.version), which does not match the tool pin")
        tool_records[key] = Provenance.as_yaml(record)
    end
    mkpath(output)
    hashfile(p) = bytes2hex(sha256(read(p)))
    mapping = Dict("run_A" => ["F3D0", "F3D1"], "run_B" => ["F3D2"])
    input_hashes = Dict{String,String}()
    for (dir, _, files) in walkdir(joinpath(root, "data", "MiSeq_SOP")), name in files
        endswith(name, ".fastq.gz") || continue
        path = joinpath(dir, name)
        input_hashes[replace(relpath(path, root), '\\' => '/')] = hashfile(path)
    end
    length(input_hashes) == 6 || error("expected exactly three paired-read samples")
    source_files = split(readchomp(`git -C $root ls-files --cached --others --exclude-standard`), '\n')
    source_hashes = Dict(replace(path, '\\' => '/') => hashfile(joinpath(root, path))
                         for path in source_files if isfile(joinpath(root, path)))
    common = Dict{String,Any}(
        "schema" => "metamanifold.pipeline-count-provenance.v1",
        "source_commit" => readchomp(`git -C $root rev-parse HEAD`),
        "source_file_hashes" => source_hashes,
        "julia" => string(VERSION), "r" => rversion, "dada2" => dada2_version,
        "embedded_r" => r_record.version,
        "r_packages" => r_record.packages,
        "tools" => tool_records,
        "manifest_sha256" => hashfile(joinpath(root, "Manifest.toml")),
        "renv_lock_sha256" => hashfile(joinpath(root, "renv.lock")),
        "tool_pins_sha256" => hashfile(joinpath(root, "config/defaults/tool_versions.yml")),
        "input_fastq_sha256" => input_hashes,
        "sample_mapping" => mapping,
        "generated_at_utc" => string(now(UTC)),
        "pipeline_started_unix" => started,
        "data_source" => "Mothur MiSeq SOP, Kozich et al. (2013); repository's three-sample subset",
        "data_use_status" => "attribution identified; independent data-rights/release approval not asserted",
        "taxonomy_claim" => "none: matrix is the pre-taxonomy DADA2 count output",
    )
    for run in sort(collect(keys(mapping)))
        project = joinpath(root, "projects", "MiSeq_SOP", run)
        matrix = joinpath(project, "dada2", "Tables", "seqtab_nochim.csv")
        _require_fresh_matrix(matrix, started)
        header = collect(first(CSV.Rows(matrix; header=false, types=String)))
        labels = String.(header[2:end])
        length(labels) == length(mapping[run]) || error("unexpected number of samples in $run")
        sample_identity = Dict(label => only(filter(id -> label == id || startswith(label, id * "_"),
                                                   mapping[run])) for label in labels)
        sort(collect(values(sample_identity))) == sort(mapping[run]) || error("sample mapping mismatch in $run")
        config = joinpath(project, "run_config.yml")
        isfile(config) || error("missing effective run configuration $config")
        directory = joinpath(output, run)
        mkpath(directory)
        cp(matrix, joinpath(directory, "counts.csv"); force=true)
        cp(config, joinpath(directory, "pipeline.yml"); force=true)
        # Keep the pipeline's per-stage probes and output hashes when available.
        attestation = joinpath(project, "attestation.yml")
        isfile(attestation) && cp(attestation, joinpath(directory, "attestation.yml"); force=true)
        provenance = copy(common)
        provenance["run"] = run
        provenance["output_sample_identity"] = sample_identity
        provenance["effective_config_sha256"] = hashfile(config)
        provenance["counts_modified_unix"] = mtime(matrix)
        provenance_path = joinpath(directory, "provenance.json")
        write(provenance_path, jcs_bytes(provenance))
        result_path = joinpath(directory, "result.json")
        biological_results(joinpath(directory, "counts.csv"), provenance_path, result_path)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: collect_pipeline_evidence.jl OUTPUT_DIRECTORY")
    collect_pipeline_evidence(ARGS[1])
end
