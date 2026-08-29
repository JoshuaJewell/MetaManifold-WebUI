module Config

# Hierarchical config cascade and per-section content-hash helpers.
#
# Config cascade (global -> study -> group(s) -> run):
#   config/defaults/pipeline.yml     <- full defaults (source of truth)
#   config/pipeline.yml              <- machine-level overrides
#   data/{study}/pipeline.yml        <- study-level overrides
#   ...intermediate dirs...
#   data/{study}/{run}/pipeline.yml  <- run-level overrides
#
# Pipeline configs live alongside the input data so that inputs are together
# in data/ and outputs remain isolated in projects/.
#
# Each override file contains only intentional changes; omitted keys are
# inherited from the nearest ancestor. At the start of each run, all levels
# are deep-merged into a single Dict and written to:
#   projects/{study}/{group}/{run}/run_config.yml
#
# Stage skip guards hash the relevant YAML section from run_config.yml so
# that any change at any level - global, study, or run - correctly
# invalidates downstream checkpoints.
#
# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export _section_stale, _write_section_hash, _begin_section, stage_run_count, _stale_keys,
       load_merged_config, write_run_config, stage_sections,
       DEFAULT_SEED, canonical_yaml_doc

    using SHA, YAML
    using ..PipelineTypes
    # Defined in Provenance, which loads first; re-exported here for callers
    # that already reach for Config when they write a config document.
    using ..Provenance: canonical_yaml_doc


    const DEFAULT_SEED = 123

    ## Pipeline stage sections
    # Maps each stage to the config sections whose change should invalidate it.
    const _STAGE_SECTIONS = Dict(
        :fastqc_multiqc => ["fastqc", "multiqc"],
        # "primers" is the sequence table write_run_config embeds from
        # primers.yml. It is the substantive input to trimming - a pair whose
        # definition changes changes the reads - so it has to invalidate this
        # stage. Without it, editing a primer sequence leaves the previously
        # trimmed reads in place and every downstream stage skips over them.
        :cutadapt => ["seed", "subsample_n", "cutadapt", "primers", "dada2.file_patterns.mode"],
        :dada2_filter_trim => ["dada2.file_patterns", "dada2.filter_trim"],
        :dada2_learn_errors => ["seed", "dada2.dada"],
        :dada2_denoise => ["seed", "dada2.dada", "dada2.merge"],
        :dada2_filter_length => ["dada2.asv"],
        :dada2_chimera_removal => ["dada2.asv", "dada2.output"],
        # "reference_database" is the databases.yml entry write_run_config embeds,
        # so a new release under the same key invalidates classification.
        :dada2_assign_taxonomy => ["seed", "dada2.taxonomy", "dada2.output", "reference_database.dada2"],
        :cdhit => ["cdhit"],
        :swarm => ["swarm"],
        :vsearch => ["vsearch", "dada2.taxonomy.database", "reference_database.vsearch"],
        :merge_taxa => ["merge_taxa", "tagging", "vsearch.enabled", "swarm.enabled", "dada2.taxonomy.enabled",
                        "reference_database.levels", "reference_database.corrections"],
    )

    function stage_sections(name::Symbol)::String
        sections = get(_STAGE_SECTIONS, name, nothing)
        isnothing(sections) && error("No pipeline stage named :$name")
        join(sections, ",")
    end

    ## Deep merge
    # Recursively merge `patch` into `base`. Dict values are merged recursively;
    # all other types (including Arrays) are replaced by the patch value.
    function _deep_merge(base::Dict, patch::Dict)
        result = copy(base)
        for (k, v) in patch
            result[k] = (haskey(result, k) && result[k] isa Dict && v isa Dict) ?
                        _deep_merge(result[k], v) : v
        end
        return result
    end

    ## Cascade path discovery
    # Return the ordered list of pipeline.yml paths participating in the cascade,
    # from least specific (defaults) to most specific (leaf project dir).
    function _cascade_paths(config_dir::String, study_dir::String, project_dir::String)
        paths = String[
            joinpath(config_dir, "defaults", "pipeline.yml"),
            joinpath(config_dir, "pipeline.yml"),
        ]

        study_norm   = normpath(abspath(study_dir))
        project_norm = normpath(abspath(project_dir))

        if study_norm == project_norm
            push!(paths, joinpath(project_norm, "pipeline.yml"))
        else
            rel   = relpath(project_norm, study_norm)
            parts = splitpath(rel)
            current = study_norm
            push!(paths, joinpath(current, "pipeline.yml"))
            for part in parts
                current = joinpath(current, part)
                push!(paths, joinpath(current, "pipeline.yml"))
            end
        end

        return paths
    end

    ## Config loading
    # Merge all config levels in `paths` (ordered global -> specific).
    # paths[1] must exist (the defaults file). Subsequent files are optional;
    # missing files and files that parse as nothing/empty are skipped.
    function load_merged_config(paths::Vector{String})
        isfile(paths[1]) || error("Default config not found: $(paths[1])")
        base = something(YAML.load_file(paths[1]), Dict())
        base isa Dict || (base = Dict())
        for path in paths[2:end]
            isfile(path) || continue
            patch = YAML.load_file(path)
            isnothing(patch) && continue
            patch isa Dict  || continue
            isempty(patch)  && continue
            base = _deep_merge(base, patch)
        end
        return base
    end

    load_merged_config(config_dir::String, study_dir::String, project_dir::String) =
        load_merged_config(_cascade_paths(config_dir, study_dir, project_dir))

    load_merged_config(project::ProjectCtx) =
        load_merged_config(project.config_dir, project.data_study_dir, project.data_dir)

    # Narrow the primer library to the pairs this run actually trims with.
    #
    # Two reasons. A run config should record only what the run used; and this block is hashed into the cutadapt stage
    # (`_STAGE_SECTIONS[:cutadapt]`), so embedding the library wholesale would
    # invalidate trimming for every project whenever any unrelated primer was
    # edited. Narrowed, only a run whose own pairs changed goes stale.
    function _run_primers(primers, merged)
        primers isa Dict || return primers
        pairs = get(primers, "Pairs", nothing)
        pairs isa AbstractVector || return primers

        wanted = get(_get_nested(merged, "cutadapt") isa Dict ?
                     _get_nested(merged, "cutadapt") : Dict(), "primer_pairs", nothing)
        wanted isa AbstractVector || return primers
        want = [string(w) for w in wanted]

        fwd_all = get(primers, "Forward", Dict())
        rev_all = get(primers, "Reverse", Dict())
        kept_pairs = Any[]
        fwd = Dict{String,Any}()
        rev = Dict{String,Any}()

        for name in want
            for entry in pairs
                entry isa Dict || continue
                haskey(entry, name) || continue
                push!(kept_pairs, entry)
                members = entry[name]
                members isa AbstractVector || continue
                # A pair names its forward first and its reverse second, but a
                # name is looked up in both tables so an unusual ordering still
                # carries its sequence into the hash.
                for m in members
                    k = string(m)
                    haskey(fwd_all, k) && (fwd[k] = fwd_all[k])
                    haskey(rev_all, k) && (rev[k] = rev_all[k])
                end
                break
            end
        end

        Dict{String,Any}("Forward" => fwd, "Reverse" => rev, "Pairs" => kept_pairs)
    end

    # The databases.yml entry the run classifies against, without the machine-local
    # paths, so its release identity (uri, checksum, levels) is part of the run config.
    function _run_reference_database(dbs_path::String, merged)
        isfile(dbs_path) || return nothing
        doc = YAML.load_file(dbs_path)
        doc isa Dict || return nothing
        dbs = get(doc, "databases", nothing)
        dbs isa Dict || return nothing
        tax = _get_nested(merged, "dada2.taxonomy")
        key = string(tax isa Dict ? get(tax, "database", "pr2") : "pr2")
        entry = get(dbs, key, nothing)
        entry isa Dict || return nothing
        out = Dict{String,Any}("key" => key)
        for (k, v) in entry
            if v isa Dict
                out[k] = Dict{String,Any}(kk => vv for (kk, vv) in v if !(kk in ("local", "remote_path")))
            else
                out[k] = v
            end
        end
        out
    end

    ## run_config.yml
    # Merge all cascade levels and write the result to {project.dir}/run_config.yml.
    # Regenerates only when a source file is newer than the existing run_config.yml.
    # Returns the path to run_config.yml.
    function write_run_config(project::ProjectCtx)
        run_config_path = joinpath(project.dir, "run_config.yml")
        primers_path    = joinpath(project.config_dir, "primers.yml")
        dbs_path        = joinpath(project.config_dir, "databases.yml")
        paths = _cascade_paths(project.config_dir, project.data_study_dir, project.data_dir)
        source_paths = filter(isfile, vcat(paths, primers_path, dbs_path))
        if !isfile(run_config_path) ||
           any(p -> isfile(p) && mtime(p) >= mtime(run_config_path), source_paths)
            merged = load_merged_config(paths)
            if isfile(primers_path)
                primers = YAML.load_file(primers_path)
                isnothing(primers) || (merged["primers"] = _run_primers(primers, merged))
            end
            ref = _run_reference_database(dbs_path, merged)
            isnothing(ref) || (merged["reference_database"] = ref)
            text = YAML.write(canonical_yaml_doc(merged))
            # This runs from pipeline stages and from the run-list HTTP route, and
            # stages re-read the file throughout a run. Writing in place truncated
            # it first, so a stage reading at that moment could load an empty or
            # partial config. Write a sibling temp file and rename it over the
            # target: the rename is atomic, so a reader sees the old file or the
            # new one, never a torn one. An unchanged config is not rewritten at
            # all, only touched so the mtime check stops re-merging it.
            if isfile(run_config_path) && read(run_config_path, String) == text
                touch(run_config_path)
            else
                tmp, io = mktemp(dirname(run_config_path); cleanup=false)
                try
                    print(io, text)
                    close(io)
                    # rename(2), not `mv(...; force=true)`: Julia's mv removes the
                    # destination first, so the path briefly does not exist and a
                    # concurrent reader gets ENOENT. rename replaces in one step.
                    Base.Filesystem.rename(tmp, run_config_path)
                catch
                    close(io)
                    rm(tmp; force=true)
                    rethrow()
                end
            end
        end
        return run_config_path
    end

    ## Section content hashes
    # Produce a stable, canonical string from a YAML-loaded value.
    # Dicts are sorted by key so insertion-order differences do not affect the hash.
    function _canonical(x)::String
        if x isa AbstractDict
            pairs_sorted = sort([(string(k), _canonical(v)) for (k, v) in x]; by = p -> p[1])
            return "{" * join(["$(p[1]):$(p[2])" for p in pairs_sorted], ",") * "}"
        elseif x isa AbstractVector
            return "[" * join(map(_canonical, x), ",") * "]"
        elseif isnothing(x)
            return "null"
        else
            return string(x)
        end
    end

    # Recursively collect all leaf (non-Dict) values into `out` as dotted-key => canonical-string pairs.
    function _collect_leaf_values!(out::Dict{String,String}, prefix::AbstractString, val)
        if val isa Dict
            for (k, v) in val
                _collect_leaf_values!(out, prefix == "" ? string(k) : "$prefix.$k", v)
            end
        else
            out[prefix] = _canonical(val)
        end
    end

    # Walk a dotted key path ("dada2.filter_trim") into a nested Dict.
    # Returns nothing if any key is missing or a non-Dict is encountered mid-path.
    function _get_nested(cfg, path::AbstractString)
        val = cfg
        for k in split(path, ".")
            val isa Dict || return nothing
            val = get(val, k, nothing)
            isnothing(val) && return nothing
        end
        return val
    end

    # Section can be a dotted path ("dada2.filter_trim") or a comma-separated
    # list of dotted paths ("dada2.dada,dada2.merge") whose canonical strings
    # are joined before hashing.
    function _section_hash(config_path::String, section::String)::String
        cfg = YAML.load_file(config_path)
        cfg isa Dict || return bytes2hex(sha256(""))
        combined = join([_canonical(_get_nested(cfg, strip(s)))
                         for s in split(section, ",")], "|")
        bytes2hex(sha256(combined))
    end

    """
        _section_stale(config_path, section, hash_file) -> Bool

    Return `true` if the named section of `config_path` has changed since
    `hash_file` was last written, or if `hash_file` does not yet exist.
    Pass `run_config.yml` as `config_path` to hash the fully-merged section.
    """
    function _section_stale(config_path::String, section::String, hash_file::String)::Bool
        !isfile(hash_file) && return true
        stored = strip(read(hash_file, String))
        return _section_hash(config_path, section) != stored
    end

    """
        _write_section_hash(config_path, section, hash_file)

    Write the current SHA-256 hash of the named section to `hash_file`.
    Call this after a stage completes successfully.
    """
    function _write_section_hash(config_path::String, section::String, hash_file::String;
                                 snapshot=nothing)
        snap = isnothing(snapshot) ? _section_snapshot(config_path, section) : snapshot
        write(hash_file, snap.hash)
        _write_values_file(hash_file, snap.values)
    end

    """
        _begin_section(config_path, section, hash_file) -> snapshot

    Called once a stage has decided to run. Removes the recorded hash, so outputs
    from an interrupted run are never taken as current, and returns the config the
    stage runs with; pass it to `_write_section_hash` as `snapshot` on success.
    """
    function _begin_section(config_path::String, section::String, hash_file::String)
        lock(_STAGE_RUNS_LOCK) do
            key = abspath(config_path)
            _STAGE_RUNS[key] = get(_STAGE_RUNS, key, 0) + 1
        end
        snap = _section_snapshot(config_path, section)
        rm(hash_file; force=true)
        rm(hash_file * ".values"; force=true)
        snap
    end

    # Stages that did work, per run config, so a caller can tell a skip from a run.
    const _STAGE_RUNS = Dict{String,Int}()
    const _STAGE_RUNS_LOCK = ReentrantLock()
    stage_run_count(config_path::String) = lock(() -> get(_STAGE_RUNS, abspath(config_path), 0), _STAGE_RUNS_LOCK)

    _section_snapshot(config_path::String, section::String) =
        (; hash=_section_hash(config_path, section), values=_section_values(config_path, section))

    """
        _stale_keys(config_path, section, hash_file) -> Vector{String}

    If the section is stale, return the individual dotted config keys that
    changed compared to the snapshot stored in `hash_file.values`.

    When the stage last ran, `_write_section_hash` stores the hash; we also
    store a companion `.values` file with the canonical per-key values. If that
    file is missing (e.g. older run), we fall back to returning the section
    names only.
    """
    function _stale_keys(config_path::String, section::String, hash_file::String)::Vector{String}
        cfg = YAML.load_file(config_path)
        cfg isa Dict || return String[]

        values_file = hash_file * ".values"
        sections = [strip(s) for s in split(section, ",")]

        # Collect current leaf keys (recursive)
        current = Dict{String,String}()
        for sec in sections
            val = _get_nested(cfg, sec)
            isnothing(val) && continue
            _collect_leaf_values!(current, sec, val)
        end

        # Read stored values from companion file
        if isfile(values_file)
            stored = Dict{String,String}()
            for line in eachline(values_file)
                idx = findfirst('=', line)
                isnothing(idx) && continue
                stored[line[1:idx-1]] = line[idx+1:end]
            end
            # Diff: keys present in current but not stored, or with different values
            changed = String[]
            for (k, v) in current
                if !haskey(stored, k) || stored[k] != v
                    push!(changed, k)
                end
            end
            # Keys removed in current config
            for k in keys(stored)
                haskey(current, k) || push!(changed, k)
            end
            return sort(changed)
        else
            # No companion file - return section names as fallback
            return sections
        end
    end

    # The leaf values of a section, recorded beside its hash so staleness can name the changed keys.
    function _section_values(config_path::String, section::String)::Dict{String,String}
        leaves = Dict{String,String}()
        cfg = YAML.load_file(config_path)
        cfg isa Dict || return leaves
        for sec in (strip(s) for s in split(section, ","))
            val = _get_nested(cfg, sec)
            isnothing(val) && continue
            _collect_leaf_values!(leaves, sec, val)
        end
        leaves
    end

    function _write_values_file(hash_file::String, leaves::Dict{String,String})
        open(hash_file * ".values", "w") do io
            for (k, v) in sort(collect(leaves); by=first)
                println(io, "$k=$v")
            end
        end
    end

end
