module ProjectSetup

export new_project

    using ..PipelineTypes

    # Walk upward from dirname(dir) to find the nearest ancestor within
    # root_project that already has fname, then check config_dir, then
    # fall back to defaults_dir.
    function _find_config_source(dir::String, fname::String,
                                  root_project::String,
                                  config_dir::String,
                                  defaults_dir::String)
        root_norm   = normpath(root_project)
        candidate   = normpath(dirname(dir))
        while candidate == root_norm ||
              startswith(candidate * "/", root_norm * "/")
            path = joinpath(candidate, fname)
            isfile(path) && return path
            candidate == root_norm && break
            candidate = normpath(dirname(candidate))
        end
        global_path = joinpath(config_dir, fname)
        isfile(global_path) && return global_path
        return joinpath(defaults_dir, fname)
    end

    """
        new_project(name; data_dir, projects_dir, config_dir)

    Bootstrap a project tree under `projects_dir/{name}/` mirroring the
    directory structure found under `data_dir/{name}/`.

    Any directory within `data_dir/{name}/` that contains `.fastq.gz` files is
    treated as a leaf run. Intermediate directories and leaf runs each receive
    copies of the stage config templates, cascading downward: each level
    inherits from the nearest ancestor that already has a given config, falling
    back to `config_dir/defaults/` when no ancestor has it.

    Config cascade (nearest-first): run -> study root -> intermediate levels
    -> `config_dir/` (global) -> `config_dir/defaults/` (ultimate fallback).
    All levels are bootstrapped at call time; `config_dir/` stage configs are
    created alongside `databases.yml` and `tools.yml`.

    Re-running `new_project` never overwrites existing configs. To reset a
    level to its parent's settings, delete that config file and re-run.

    Returns a `Vector{ProjectCtx}`, one per leaf run.
    """
    function new_project(name::String;
                         data_dir::String     = "./data",
                         projects_dir::String = "./projects",
                         config_dir::String   = "./config")

        defaults_dir = joinpath(config_dir, "defaults")
        root_data    = joinpath(data_dir, name)
        root_project = joinpath(projects_dir, name)

        isdir(root_data) || error("Data directory not found: $root_data")

        stage_configs = ("dada2.yml", "cutadapt.yml", "vsearch.yml",
                         "cdhit.yml", "merge_taxa.yml")

        # Bootstrap global configs from defaults if missing.
        for fname in ("databases.yml", "tools.yml", stage_configs...)
            dst = joinpath(config_dir, fname)
            if !isfile(dst)
                src = joinpath(defaults_dir, fname)
                isfile(src) || error("Default template not found: $src")
                cp(src, dst)
                @info "Created: $dst"
            end
        end

        # Find leaf directories (those containing .fastq.gz), relative to root_data.
        leaf_relpaths = String[]
        for (dirpath, _, files) in walkdir(root_data)
            if any(f -> endswith(f, ".fastq.gz"), files)
                push!(leaf_relpaths, relpath(dirpath, root_data))
            end
        end
        isempty(leaf_relpaths) && error("No .fastq.gz files found under $root_data")

        # Collect all directories to bootstrap (study root, intermediates, leaves).
        # root_project is always included so it sits in the cascade between
        # config/ and any deeper levels.
        all_dirs = Set{String}([root_project])
        for rp in leaf_relpaths
            rp == "." && continue
            parts = splitpath(rp)
            for i in 1:length(parts)
                push!(all_dirs, joinpath(root_project, parts[1:i]...))
            end
        end

        # Process top-down so each level can serve as source for its children.
        sorted_dirs = sort(collect(all_dirs), by = d -> length(splitpath(d)))

        for dir in sorted_dirs
            mkpath(dir)
            for fname in stage_configs
                dst = joinpath(dir, fname)
                if isfile(dst)
                    @info "Exists (skipping): $dst"
                    continue
                end
                src = _find_config_source(dir, fname, root_project, config_dir, defaults_dir)
                isfile(src) || error("Default template not found: $src")
                cp(src, dst)
                @info "Created: $dst"
            end
        end

        # Return one ProjectCtx per leaf.
        projects = ProjectCtx[]
        for rp in leaf_relpaths
            proj_dir  = rp == "." ? root_project : joinpath(root_project, rp)
            fastq_dir = rp == "." ? root_data    : joinpath(root_data,    rp)
            push!(projects, ProjectCtx(proj_dir, config_dir, fastq_dir))
        end
        @info "$(length(projects)) project(s) ready under $root_project"
        return projects
    end

end
