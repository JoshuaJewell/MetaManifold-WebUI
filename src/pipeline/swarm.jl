module OTUPipeline

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export swarm

    using CSV, DataFrames, YAML, Logging
    using ..PipelineTypes
    using ..PipelineGraph
    using ..PipelineLog
    using ..Config
    using ..Tools: tool_bin

    function _run_logged(cmd_str::String, log_path::String; mode::String="w")
        try
            open(log_path, mode) do io
                run(pipeline(`bash -lc $cmd_str`; stdout=io, stderr=io))
            end
        catch e
            isfile(log_path) && print(stderr, read(log_path, String))
            rethrow()
        end
    end

    function _swarm_args(cfg::Dict)::String
        parts = String[]
        push!(parts, "-d $(get(cfg, "differences", 1))")
        threads = get(cfg, "threads", 0)
        threads > 0 && push!(parts, "-t $threads")
        extra = strip(get(cfg, "optional_args", ""))
        isempty(extra) || push!(parts, extra)
        join(parts, " ")
    end

    # Latest mtime of trimmed FASTQ files.
    function _trimmed_mtime(trimmed::TrimmedReads)
        files = filter(f -> endswith(f, ".fastq.gz"), readdir(trimmed.dir))
        isempty(files) ? 0.0 : maximum(mtime(joinpath(trimmed.dir, f)) for f in files)
    end

    ## Step 1: merge paired-end reads per sample into per-sample FASTA files.
    # Returns the list of sample stems that were found.
    function _merge_pairs!(merged_dir::String, trimmed::TrimmedReads, cfg::Dict,
                           vsearch_bin::String, log_dir::String)
        mkpath(merged_dir)
        mode     = get(get(cfg, "file_patterns", Dict()), "mode", "paired")
        minovlen = get(cfg, "fastq_minovlen", 20)
        r1_files = filter(f -> occursin(r"_R1_trimmed\.fastq\.gz$", f), readdir(trimmed.dir))
        samples  = [replace(f, r"_R1_trimmed\.fastq\.gz$" => "") for f in r1_files]
        for sample in samples
            r1        = joinpath(trimmed.dir, sample * "_R1_trimmed.fastq.gz")
            out_fasta = joinpath(merged_dir, sample * ".fasta")
            isfile(out_fasta) && mtime(out_fasta) > mtime(r1) && continue
            log_path = joinpath(log_dir, "merge_$(sample).log")
            if mode == "paired"
                r2  = joinpath(trimmed.dir, sample * "_R2_trimmed.fastq.gz")
                cmd = "$(vsearch_bin) --fastq_mergepairs $(r1) --reverse $(r2) --fastaout $(out_fasta) --fastq_minovlen $(minovlen) --quiet"
            else
                cmd = "$(vsearch_bin) --fastq_filter $(r1) --fastaout $(out_fasta) --quiet"
            end
            _run_logged(cmd, log_path)
        end
        return samples
    end

    ## Step 2: pool all per-sample FASTAs, dereplicate, sort by abundance.
    function _derep_sort!(swarm_dir::String, merged_dir::String, samples::Vector{String},
                          min_abundance::Int, vsearch_bin::String, log_dir::String)
        pooled = joinpath(swarm_dir, "pooled.fasta")
        derep  = joinpath(swarm_dir, "derep.fasta")
        sorted = joinpath(swarm_dir, "sorted.fasta")
        open(pooled, "w") do io
            for sample in samples
                fa = joinpath(merged_dir, sample * ".fasta")
                isfile(fa) && write(io, read(fa, String))
            end
        end
        _run_logged("$(vsearch_bin) --derep_fulllength $(pooled) --output $(derep) --sizeout --quiet",
                    joinpath(log_dir, "derep.log"))
        _run_logged("$(vsearch_bin) --sortbysize $(derep) --output $(sorted) --minsize $(min_abundance) --quiet",
                    joinpath(log_dir, "sort.log"))
        return sorted
    end

    ## Step 3: chimera filtering with vsearch uchime_denovo (optional).
    function _chimera_filter!(sorted::String, swarm_dir::String, vsearch_bin::String, log_dir::String)
        nochim = joinpath(swarm_dir, "nochim.fasta")
        _run_logged("$(vsearch_bin) --uchime_denovo $(sorted) --nonchimeras $(nochim) --quiet",
                    joinpath(log_dir, "chimera.log"))
        return nochim
    end

    ## Step 3b: remove sequences containing ambiguous bases (N).
    # Swarm rejects any IUPAC ambiguity code; vsearch --fastx_filter drops them cleanly.
    function _filter_ns!(input::String, swarm_dir::String, vsearch_bin::String, log_dir::String)
        noN = joinpath(swarm_dir, "noN.fasta")
        _run_logged("$(vsearch_bin) --fastx_filter $(input) --fastaout $(noN) --fastq_maxns 0 --sizeout --quiet",
                    joinpath(log_dir, "filter_ns.log"))
        return noN
    end

    ## Step 4: SWARM clustering. Seeds are renamed otu1, otu2, ... for clean IDs.
    function _cluster!(input_fasta::String, swarm_dir::String, args::String,
                       swarm_bin::String, log_dir::String)
        seeds_raw = joinpath(swarm_dir, "seeds_raw.fasta")
        seeds     = joinpath(swarm_dir, "seeds.fasta")
        swarm_txt = joinpath(swarm_dir, "swarm.txt")
        _run_logged("$(swarm_bin) -z $(args) $(input_fasta) --seeds $(seeds_raw) --uclust-file $(swarm_txt)",
                    joinpath(log_dir, "swarm.log"))
        # Rename to otu1, otu2, ... for consistent IDs across count table and taxonomy search.
        _run_logged("awk 'BEGIN{n=0}/^>/{printf \">otu%d\\n\",++n;next}{print}' $(seeds_raw) > $(seeds)",
                    joinpath(log_dir, "rename.log"))
        return seeds
    end

    ## Step 5: map per-sample merged reads back to OTU seeds, build count table CSV.
    function _build_count_table!(merged_dir::String, seeds::String, samples::Vector{String},
                                  identity::Float64, vsearch_bin::String,
                                  swarm_dir::String, log_dir::String)
        counts_dir = joinpath(swarm_dir, "counts")
        mkpath(counts_dir)
        otu_ids = String[]
        open(seeds, "r") do io
            for line in eachline(io)
                startswith(line, ">") && push!(otu_ids, strip(line[2:end]))
            end
        end
        sample_counts = Dict{String, Dict{String,Int}}()
        for sample in samples
            fa = joinpath(merged_dir, sample * ".fasta")
            isfile(fa) || continue
            hits_file = joinpath(counts_dir, sample * "_hits.txt")
            cmd = "$(vsearch_bin) --usearch_global $(fa) --db $(seeds) --userout $(hits_file) --userfields target --maxaccepts 1 --id $(identity) --quiet"
            _run_logged(cmd, joinpath(log_dir, "count_$(sample).log"))
            counts = Dict{String,Int}()
            open(hits_file, "r") do io
                for line in eachline(io)
                    otu = strip(line)
                    isempty(otu) || (counts[otu] = get(counts, otu, 0) + 1)
                end
            end
            sample_counts[sample] = counts
        end
        df = DataFrame(SeqName = otu_ids)
        for sample in samples
            sc = get(sample_counts, sample, Dict{String,Int}())
            df[!, sample] = [get(sc, otu, 0) for otu in otu_ids]
        end
        otu_table = joinpath(swarm_dir, "otu_table.csv")
        CSV.write(otu_table, df)
        @info "SWARM: $(length(otu_ids)) OTUs written to $otu_table"
        return otu_table
    end

    """
        swarm(project, trimmed; swarm_bin, vsearch_bin) -> Union{OTUResult, Nothing}

    Run the SWARM OTU clustering pipeline on trimmed reads:
      merge pairs -> pool -> dereplicate -> sort -> [chimera filter] -> cluster -> count table.

    Outputs are written to `{project.dir}/swarm/`. Returns an OTUResult pointing to
    `seeds.fasta` (OTU representatives) and `otu_table.csv` (per-sample counts).

    Skips the entire pipeline when outputs are up to date and config is unchanged.
    Per-sample merge steps are individually skipped when already current.
    """
    function swarm(project::ProjectCtx, trimmed::TrimmedReads;
                   swarm_bin   = tool_bin("swarm"),
                   vsearch_bin = tool_bin("vsearch"))
        if isnothing(Sys.which(swarm_bin))
            @warn "swarm binary not found ('$swarm_bin') - skipping OTU pipeline. Install swarm or set path in config/tools.yml."
            return nothing
        end
        config_path = write_run_config(project)
        cfg         = get(YAML.load_file(config_path), "swarm", Dict())
        swarm_dir   = joinpath(project.dir, "swarm")
        merged_dir  = joinpath(swarm_dir, "merged")
        log_dir     = joinpath(swarm_dir, "logs")
        hash_file   = joinpath(swarm_dir, "config.hash")
        seeds       = joinpath(swarm_dir, "seeds.fasta")
        otu_table   = joinpath(swarm_dir, "otu_table.csv")

        t_mtime = _trimmed_mtime(trimmed)
        if isfile(seeds) && isfile(otu_table) &&
           !_section_stale(config_path, stage_sections(:swarm), hash_file) &&
           mtime(seeds) > t_mtime && mtime(otu_table) > t_mtime
            @info "SWARM: skipping - outputs up to date in $swarm_dir"
            return OTUResult(seeds, otu_table)
        end

        mkpath(merged_dir); mkpath(log_dir)

        args          = _swarm_args(cfg)
        min_abundance = get(cfg, "min_abundance", 2)
        do_chimera    = get(cfg, "chimera_check", true)
        identity      = Float64(get(cfg, "identity", 0.97))

        @info "SWARM: merging paired reads"
        samples = _merge_pairs!(merged_dir, trimmed, cfg, vsearch_bin, log_dir)

        @info "SWARM: dereplicating and sorting (minsize=$min_abundance)"
        sorted = _derep_sort!(swarm_dir, merged_dir, samples, min_abundance, vsearch_bin, log_dir)

        after_chimera = if do_chimera
            @info "SWARM: chimera filtering"
            _chimera_filter!(sorted, swarm_dir, vsearch_bin, log_dir)
        else
            sorted
        end

        @info "SWARM: removing ambiguous-base sequences"
        cluster_input = _filter_ns!(after_chimera, swarm_dir, vsearch_bin, log_dir)

        @info "SWARM: clustering"
        seeds = _cluster!(cluster_input, swarm_dir, args, swarm_bin, log_dir)

        @info "SWARM: building count table (identity=$identity)"
        otu_table = _build_count_table!(merged_dir, seeds, samples, identity,
                                         vsearch_bin, swarm_dir, log_dir)

        _write_section_hash(config_path, stage_sections(:swarm), hash_file)
        pipeline_log(project, "SWARM complete")
        log_written(project, seeds)
        log_written(project, otu_table)
        return OTUResult(seeds, otu_table)
    end

end # module OTUPipeline
