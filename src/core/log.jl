module PipelineLog

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export pipeline_log, log_written, reset_log, finalise_log, write_combined_log

    using SHA, Dates
    using ..PipelineTypes

    function reset_log(dir::String)
        path = joinpath(dir, "pipeline.log")
        open(path, "w") do io
            println(io, "Pipeline log started ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        end
    end
    reset_log(project::ProjectCtx) = reset_log(project.dir)

    function pipeline_log(dir::String, msg::String)
        open(joinpath(dir, "pipeline.log"), "a") do io
            println(io, Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "  ", msg)
        end
    end
    pipeline_log(project::ProjectCtx, msg::String) = pipeline_log(project.dir, msg)

    function log_written(dir::String, path::String)
        h = open(io -> bytes2hex(sha256(io)), path, "r")
        rel = "./" * relpath(path, dir)
        pipeline_log(dir, "Written: $rel  sha256:$h")
    end
    log_written(project::ProjectCtx, path::String) = log_written(project.dir, path)

    # Tool-level log/stats/hash files to embed, relative to the project dir.
    # Each entry is (stage_label, relative_path).  Only existing files are included.
    const _TOOL_LOG_FILES = [
        ("QC",       "QC/logs/fastqc.log"),
        ("QC",       "QC/logs/multiqc.log"),
        ("cutadapt", "cutadapt/logs/cutadapt_primer_trimming_stats.txt"),
        ("cutadapt", "cutadapt/logs/cutadapt_trimmed_percentage.txt"),
        ("cutadapt", "cutadapt/config.hash"),
        ("DADA2",    "dada2/Logs/prefilter_qc.log"),
        ("DADA2",    "dada2/Logs/filter_trim.log"),
        ("DADA2",    "dada2/Logs/learn_errors.log"),
        ("DADA2",    "dada2/Logs/denoise.log"),
        ("DADA2",    "dada2/Logs/filter_length.log"),
        ("DADA2",    "dada2/Logs/chimera_removal.log"),
        ("DADA2",    "dada2/Logs/assign_taxonomy.log"),
        ("DADA2",    "dada2/Checkpoints/filter_trim.hash"),
        ("DADA2",    "dada2/Checkpoints/learn_errors.hash"),
        ("DADA2",    "dada2/Checkpoints/denoise.hash"),
        ("DADA2",    "dada2/Checkpoints/filter_length.hash"),
        ("DADA2",    "dada2/Checkpoints/chimera_removal.hash"),
        ("DADA2",    "dada2/Checkpoints/assign_taxonomy.hash"),
        ("VSEARCH",  "vsearch/logs/vsearch.log"),
        ("VSEARCH",  "vsearch/config.hash"),
        ("cd-hit",   "cdhit/logs/cdhit.log"),
        ("cd-hit",   "cdhit/config.hash"),
        ("merge",    "merged/config.hash"),
    ]

    """
        finalise_log(project::ProjectCtx)

    Append the contents of all tool-level log, stats, and hash files into
    the run's `pipeline.log`.  Call once per run after all stages complete.
    """
    function finalise_log(project::ProjectCtx)
        pipeline_log(project, "--- Tool logs appended below ---")
        prev_stage = ""
        for (stage, relpath_) in _TOOL_LOG_FILES
            abspath_ = joinpath(project.dir, relpath_)
            isfile(abspath_) || continue
            if stage != prev_stage
                pipeline_log(project, "")  # blank separator
                prev_stage = stage
            end
            _append_tool_file(project.dir, relpath_, abspath_)
        end
        pipeline_log(project, "--- End of tool logs ---")
    end

    # Append a single tool file into the pipeline.log with a header/footer.
    function _append_tool_file(dir::String, relpath_::String, abspath_::String)
        contents = read(abspath_, String)
        open(joinpath(dir, "pipeline.log"), "a") do io
            println(io, ">>>>>> ", relpath_)
            print(io, contents)
            if !isempty(contents) && !endswith(contents, '\n')
                println(io)
            end
            println(io, "<<<<<< ", relpath_)
        end
    end

    """
        write_combined_log(projects; study_dir)

    Finalise each run's log (appending tool-level files), then write a
    combined log to `{study_dir}/combined_pipeline.log` that merges all
    per-run `pipeline.log` files, each run's `run_config.yml`, plus any
    group- and study-level logs.
    """
    function write_combined_log(projects::Vector{ProjectCtx};
                                study_dir::String = projects[1].study_dir)
        # Finalise each run's pipeline.log with tool-level contents.
        for project in projects
            finalise_log(project)
        end

        out_path = joinpath(study_dir, "combined_pipeline.log")
        sep = "=" ^ 72

        open(out_path, "w") do io
            println(io, sep)
            println(io, "  Combined pipeline log")
            println(io, "  Generated: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
            println(io, "  Study: ", study_dir)
            println(io, "  Runs:  ", length(projects))
            println(io, sep)

            # Per-run sections
            for project in projects
                run_name = relpath(project.dir, study_dir)
                println(io, "\n", sep)
                println(io, "  RUN: ", run_name)
                println(io, sep)

                # run_config.yml
                config_path = joinpath(project.dir, "run_config.yml")
                if isfile(config_path)
                    println(io, "\n--- run_config.yml ---")
                    print(io, read(config_path, String))
                    println(io)
                end

                # pipeline.log (now includes tool logs)
                log_path = joinpath(project.dir, "pipeline.log")
                if isfile(log_path)
                    println(io, "--- pipeline.log ---")
                    print(io, read(log_path, String))
                end
            end

            # Group-level logs (directories between runs and study root)
            group_dirs = Set{String}()
            for project in projects
                gdir = dirname(project.dir)
                while gdir != study_dir && startswith(gdir, study_dir)
                    push!(group_dirs, gdir)
                    gdir = dirname(gdir)
                end
            end

            for gdir in sort(collect(group_dirs))
                log_path = joinpath(gdir, "pipeline.log")
                isfile(log_path) || continue
                gname = relpath(gdir, study_dir)
                println(io, "\n", sep)
                println(io, "  GROUP: ", gname)
                println(io, sep)
                println(io, "\n--- pipeline.log ---")
                print(io, read(log_path, String))
            end

            # Study-level log
            study_log = joinpath(study_dir, "pipeline.log")
            if isfile(study_log)
                println(io, "\n", sep)
                println(io, "  STUDY")
                println(io, sep)
                println(io, "\n--- pipeline.log ---")
                print(io, read(study_log, String))
            end
        end

        @info "Written: $out_path"
        return out_path
    end
end
