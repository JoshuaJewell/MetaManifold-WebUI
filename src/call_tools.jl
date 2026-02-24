module Tools

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export cutadapt, vsearch, multiqc, cdhit

    using YAML
    using Logging
    using ..PipelineTypes
    using ..PipelineLog
    using ..Config

    # Run cmd_str via bash, capturing stdout+stderr to log_path.
    # On failure, prints the log to stderr before rethrowing so the error is visible.
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

    ## Tools loading from config
    """
        load_tools(config_path) -> Dict{String,String}

    Read config/tools.yml and return a Dict of tool name => resolved path.
    If the file does not exist, returns an empty Dict so tools fall back to PATH.
    Paths with `@` are SSH remote paths (user@host:/path), the calling module is
    responsible for routing those calls via SSH.
    """
    function load_tools(config_path = joinpath(@__DIR__, "..", "config", "tools.yml"))
        isfile(config_path) || return Dict{String,String}()
        data = YAML.load_file(config_path)
        tools = Dict{String,String}()
        for (name, info) in data
            path = info isa Dict ? get(info, "path", nothing) : nothing
            path !== nothing && (tools[name] = string(path))
        end
        tools
    end

    tool_bin(key) = get(_tools, key, key)

    const _tools = load_tools()

    ## cutadapt
    function get_primers(input, primers_path, seen_fwd::Vector{String}, seen_rev::Vector{String})
        data = YAML.load_file(primers_path)

        for pair in data["Pairs"]
            if haskey(pair, input)
                pair_names = pair[input]

                forward = data["Forward"][pair_names[1]]
                reverse = data["Reverse"][pair_names[2]]

                output_forward = in(forward, seen_fwd) ? "x" : forward
                output_reverse = in(reverse, seen_rev) ? "x" : reverse

                if output_forward != "x"
                    push!(seen_fwd, output_forward)
                end
                if output_reverse != "x"
                    push!(seen_rev, output_reverse)
                end

                return (output_forward, output_reverse)
            end
        end

        path = pwd() * "/" * primers_path
        return "Invalid input. Please specify a valid primer pair in '$path'."
    end

    # Formats cutadapt's arguments for the primer sets requested
    function get_primer_args(primer_pairs, primers_path)
        seen_fwd = String[]
        seen_rev = String[]
        args = ""

        for pair in primer_pairs
            pair_tuple = get_primers(pair, primers_path, seen_fwd, seen_rev)
            forward_primer = pair_tuple[1]
            reverse_primer = pair_tuple[2]

            args *= join([flag * primer * " " for (flag, primer) in (("-g ", forward_primer), ("-G ", reverse_primer)) if primer != "x"], "")
        end

        return chop(args)
    end

    function run_cutadapt(primer_args, optional_args, fastq_in_dir, cutadapt_dir, cutadapt_bin)
        samples = []

        log_dir          = joinpath(cutadapt_dir, "logs")
        stats_path       = joinpath(log_dir, "cutadapt_primer_trimming_stats.txt")
        summary_path     = joinpath(log_dir, "cutadapt_trimmed_percentage.txt")
        stats_basename   = basename(stats_path)
        summary_basename = basename(summary_path)

        isdir(log_dir) || mkpath(log_dir)

        # Filter .fastq.gz files
        for f in filter(x->occursin(r"_R1.*fastq\.gz$", x), readdir(fastq_in_dir))
            push!(samples, split(f, '_')[1])
        end

        @info("cutadapt running with arguments: $primer_args $optional_args.")

        nsamples = length(samples)
        for (i, sample) in enumerate(samples)
            inputR1 = joinpath(fastq_in_dir, sample * "_*_L001_R1_001.fastq.gz")
            inputR2 = joinpath(fastq_in_dir, sample * "_*_L001_R2_001.fastq.gz")

            outputR1 = joinpath(cutadapt_dir, sample * "_R1_trimmed.fastq.gz")
            outputR2 = joinpath(cutadapt_dir, sample * "_R2_trimmed.fastq.gz")

            @info("On sample $i/$nsamples ($sample).")
            cutadapt_cmd = "$cutadapt_bin $primer_args $optional_args -o $outputR1 -p $outputR2 $inputR1 $inputR2"

            open(stats_path, "a") do io
                run(pipeline(`bash -lc $cutadapt_cmd`; stdout=io, stderr=io))
            end
        end

        samples_str = join(samples, " ")
        cmd = "paste <(printf \"%s\\n\" $samples_str) " *
            "<(grep \"passing\" $stats_basename | cut -f3 -d \"(\" | tr -d \")\") " *
            "<(grep \"filtered\" $stats_basename | cut -f3 -d \"(\" | tr -d \")\") " *
            "> $summary_basename"

        cd(log_dir) do
            run(pipeline(`bash -lc $cmd`))
        end

        @info("cutadapt complete. Output available in $cutadapt_dir.")
        return cutadapt_dir
    end

    """
        cutadapt(primer_pairs, primers_path, fastq_in_dir, cutadapt_dir; optional_args, cutadapt_bin)

    Requires `cutadapt` installed. Runs `cutadapt` with as many primer pairs as necessary
    (useful for multiplex) and any optional parameters. Primer arguments are determined from YAML file in format:
    ``` YAML
        Forward:
            PrimerF: "CCAGCASCYGCGGTAATTCC"

        Reverse:
            Primer1R: "ACTTTCGTTCTTGATYRA"
            Primer2R: "DCTKTCGTYCTTGATYRA"

        Pairs:
            - PrimerPair1:
                - PrimerF
                - Primer1R
            - PrimerPair2:
                - PrimerF
                - Primer2R
    ```
    Where a primer pair contains the same forward or reverse as another specified pair, these are deduplicated by name rather than sequence. For example:
    ```julia
    primer_pairs = ["PrimerPair1", "PrimerPair2"]
    ```
    This would execute `cutadapt` with arguments `-g CCAGCASCYGCGGTAATTCC -G ACTTTCGTTCTTGATYRA -G DCTKTCGTYCTTGATYRA`.

    ## Arguments
    - `primer_pairs`: Specify primers in string array format.
    - `primers_path`: Specify path of YAML file for primers.
    - `fastq_in_dir`: Specify path of *_*_L001_R1_001.fastq.gz and *_*_L001_R2_001.fastq.gz files for cutadapt.
    - `cutadapt_dir`: Specify an output directory for trimmed .fastq.gz and logs.

    ## Keyword Arguments
    - `optional_args` (optional, default: "-m 200 --discard-untrimmed"): Specify additional arguments passed to `cutadapt` command.
    - `cutadapt_bin` (optional, default: "cutadapt"): Path to the cutadapt binary. Defaults to PATH lookup. Set via `config/tools.yml` and `load_tools()` in main.jl.

    """
    function cutadapt(
        primer_pairs,
        primers_path,
        fastq_in_dir,
        cutadapt_dir;
        optional_args = "-m 200 --discard-untrimmed",
        cutadapt_bin  = tool_bin("cutadapt")
    )
        if isdir(cutadapt_dir)
            trimmed = filter(f -> endswith(f, "_trimmed.fastq.gz"), readdir(cutadapt_dir))
            if !isempty(trimmed) &&
               all(f -> mtime(joinpath(cutadapt_dir, f)) > mtime(primers_path), trimmed)
                @info "Skipping cutadapt: trimmed reads up to date in $cutadapt_dir"
                return TrimmedReads(cutadapt_dir)
            end
        end
        run_cutadapt(
            get_primer_args(primer_pairs, primers_path),
            optional_args,
            fastq_in_dir,
            cutadapt_dir,
            cutadapt_bin
        )
        return TrimmedReads(cutadapt_dir)
    end

    function cutadapt(project::ProjectCtx;
                      optional_args::Union{String,Nothing} = nothing,
                      cutadapt_bin = tool_bin("cutadapt"))
        config_path   = write_run_config(project)
        primers_path  = joinpath(project.config_dir, "primers.yml")
        cfg           = get(YAML.load_file(config_path), "cutadapt", Dict())
        primer_pairs  = cfg["primer_pairs"]
        optional_args = isnothing(optional_args) ?
                        get(cfg, "optional_args", "-m 200 --discard-untrimmed") : optional_args
        cutadapt_dir  = joinpath(project.dir, "cutadapt")
        hash_file     = joinpath(cutadapt_dir, "config.hash")
        if isdir(cutadapt_dir)
            trimmed = filter(f -> endswith(f, "_trimmed.fastq.gz"), readdir(cutadapt_dir))
            if !isempty(trimmed) &&
               !_section_stale(config_path, "cutadapt", hash_file) &&
               all(f -> mtime(joinpath(cutadapt_dir, f)) > mtime(primers_path), trimmed) &&
               all(f -> filesize(joinpath(cutadapt_dir, f)) > 20, trimmed)  # skip empty gzips
                @info "Skipping cutadapt: trimmed reads up to date in $cutadapt_dir"
                return TrimmedReads(cutadapt_dir)
            end
        end
        mkpath(cutadapt_dir)
        run_cutadapt(get_primer_args(primer_pairs, primers_path), optional_args,
                     project.data_dir, cutadapt_dir, cutadapt_bin)
        _write_section_hash(config_path, "cutadapt", hash_file)
        pipeline_log(project, "cutadapt complete")
        return TrimmedReads(cutadapt_dir)
    end

    ## FastQC / MultiQC
    """
        multiqc(fastq_in_dir, qc_dir; fastqc_args, multiqc_args, fastqc_bin, multiqc_bin)

    Run FastQC on all FASTQ files in `fastq_in_dir`, saving individual reports to
    `qc_dir/fastqc/`, then aggregate them into a MultiQC summary report at `qc_dir/`.

    ## Arguments
    - `fastq_in_dir`: Directory containing `.fastq` or `.fastq.gz` files.
    - `qc_dir`: Parent output directory. FastQC reports go to `qc_dir/fastqc/`;
      the MultiQC summary report goes to `qc_dir/`.

    ## Keyword Arguments
    - `fastqc_args` (optional, default: `"-t 20 --extract --delete"`): Additional arguments passed to `fastqc`.
    - `multiqc_args` (optional, default: `""`): Additional arguments passed to `multiqc`.
    - `fastqc_bin` (optional): Path to the fastqc binary. Set via `config/tools.yml`.
    - `multiqc_bin` (optional): Path to the multiqc binary. Set via `config/tools.yml`.
    """
    function multiqc(fastq_in_dir, qc_dir;
                     fastqc_args  = "-t 20 --extract --delete",
                     multiqc_args = "",
                     fastqc_bin   = tool_bin("fastqc"),
                     multiqc_bin  = tool_bin("multiqc"))
        report = joinpath(qc_dir, "multiqc_report.html")
        if isfile(report)
            fastqs = filter(f -> occursin(r"\.fastq(\.gz)?$", f), readdir(fastq_in_dir))
            if !isempty(fastqs) && all(f -> mtime(report) > mtime(joinpath(fastq_in_dir, f)), fastqs)
                @info "Skipping multiqc: $report up to date"
                return
            end
        end

        fastqc_dir = joinpath(qc_dir, "fastqc")
        log_dir    = joinpath(qc_dir, "logs")
        mkpath(fastqc_dir)
        mkpath(log_dir)

        fastqc_log  = joinpath(log_dir, "fastqc.log")
        multiqc_log = joinpath(log_dir, "multiqc.log")

        @info "FastQC running on $fastq_in_dir"
        fqc_cmd = "$fastqc_bin $fastq_in_dir/*.fastq* -o $fastqc_dir $fastqc_args"
        _run_logged(fqc_cmd, fastqc_log)
        @info "FastQC complete. Output: $fastqc_dir  Log: $fastqc_log"

        @info "MultiQC running on $fastqc_dir"
        mqc_cmd = "$multiqc_bin $fastqc_dir -o $qc_dir $multiqc_args"
        _run_logged(mqc_cmd, multiqc_log)
        @info "MultiQC complete. Output: $qc_dir  Log: $multiqc_log"
    end

    ## vsearch
    """
        vsearch(fasta_in_dir, reference_database; optional_args = "--id 0.75 --query_cov 0.8", vsearch_bin = "vsearch")

    Requires `vsearch` installed. Runs `vsearch` command with any optional parameters to perform taxonomy assignment by local alignment against specified database.

    ## Arguments
    - `fasta_in_dir`: Specify path of fasta files output by DADA2 pipeline.
    - `reference_database` (default: "./databases"): Specify path of reference database.

    ## Keyword Arguments
    - `optional_args` (optional, default: "--id 0.75 --query_cov 0.8"): Specify additional arguments passed to `vsearch` command.
    - `vsearch_bin` (optional, default: "vsearch"): Path to the vsearch binary. Defaults to PATH lookup. Set via `config/tools.yml` and `load_tools()` in main.jl.

    """
    function vsearch(fasta_in_dir, reference_database, vsearch_dir; optional_args = "--id 0.75 --query_cov 0.8", vsearch_bin = tool_bin("vsearch"))
        log_dir  = joinpath(vsearch_dir, "logs")
        mkpath(vsearch_dir)
        mkpath(log_dir)
        outfile  = joinpath(vsearch_dir, "taxonomy.tsv")
        log_path = joinpath(log_dir, "vsearch.log")
        @info "VSEARCH running: $fasta_in_dir against $(basename(reference_database))"
        # --userout with query+target+id captures the full FASTA header (including
        # description after the space), unlike --blast6out which truncates at the first space.
        # This is required for databases like SILVA where taxonomy is in the description field.
        cmd = "$vsearch_bin --usearch_global $fasta_in_dir --db $reference_database --userout $outfile --userfields query+target+id $optional_args"
        _run_logged(cmd, log_path)
        @info "VSEARCH complete. Output: $outfile  Log: $log_path"
    end

    function vsearch(input::HasFasta, reference_database::String, vsearch_dir::String;
                     optional_args = "--id 0.75 --query_cov 0.8",
                     vsearch_bin   = tool_bin("vsearch"))
        tsv = joinpath(vsearch_dir, "taxonomy.tsv")
        if isfile(tsv) && mtime(tsv) > mtime(input.fasta)
            @info "Skipping vsearch: $tsv up to date"
            return TaxonomyHits(tsv)
        end
        vsearch(input.fasta, reference_database, vsearch_dir;
                optional_args, vsearch_bin)
        return TaxonomyHits(tsv)
    end

    function vsearch(project::ProjectCtx, input::HasFasta, reference_database::String;
                     optional_args::Union{String,Nothing} = nothing,
                     vsearch_bin = tool_bin("vsearch"))
        config_path   = write_run_config(project)
        cfg           = get(YAML.load_file(config_path), "vsearch", Dict())
        optional_args = isnothing(optional_args) ?
                        get(cfg, "optional_args", "--id 0.75 --query_cov 0.8") : optional_args
        vsearch_dir   = joinpath(project.dir, "vsearch")
        tsv           = joinpath(vsearch_dir, "taxonomy.tsv")
        hash_file     = joinpath(vsearch_dir, "config.hash")
        if isfile(tsv) &&
           !_section_stale(config_path, "vsearch", hash_file) &&
           mtime(tsv) > mtime(input.fasta)
            @info "Skipping vsearch: $tsv up to date"
            return TaxonomyHits(tsv)
        end
        vsearch(input.fasta, reference_database, vsearch_dir; optional_args, vsearch_bin)
        _write_section_hash(config_path, "vsearch", hash_file)
        pipeline_log(project, "VSEARCH: $(input.fasta) against $(basename(reference_database))")
        log_written(project, tsv)
        return TaxonomyHits(tsv)
    end

    ## cd-hit-est
    """
        cdhit(fasta_in, cdhit_dir; optional_args = "-c 0.9", cdhit_bin)

    Run cd-hit-est on a FASTA file, writing the clustered output to `cdhit_dir`.
    Returns the path to the output FASTA.

    ## Arguments
    - `fasta_in`: Path to the input FASTA file.
    - `cdhit_dir`: Output directory for the clustered FASTA and `.clstr` file.

    ## Keyword Arguments
    - `optional_args` (optional, default: `"-c 0.9"`): Additional arguments passed to `cd-hit-est`.
    - `cdhit_bin` (optional): Path to the cd-hit-est binary. Set via `config/tools.yml`.
    """
    function cdhit(fasta_in, cdhit_dir; optional_args = "-c 0.9", cdhit_bin = tool_bin("cd_hit_est"))
        log_dir  = joinpath(cdhit_dir, "logs")
        mkpath(cdhit_dir)
        mkpath(log_dir)
        fasta_out = joinpath(cdhit_dir, basename(fasta_in))
        log_path  = joinpath(log_dir, "cdhit.log")
        @info "cd-hit-est running on $fasta_in"
        cmd = "$cdhit_bin -i $fasta_in -o $fasta_out $optional_args"
        _run_logged(cmd, log_path)
        @info "cd-hit-est complete. Output: $fasta_out  Log: $log_path"
        return fasta_out
    end

    function cdhit(input::ASVResult, cdhit_dir::String;
                   optional_args = "-c 0.9",
                   cdhit_bin     = tool_bin("cd_hit_est"))
        new_fasta = joinpath(cdhit_dir, basename(input.fasta))
        if isfile(new_fasta) && mtime(new_fasta) > mtime(input.fasta)
            @info "Skipping cdhit: $new_fasta up to date"
            return ASVResult(new_fasta, input.count_table, input.taxonomy)
        end
        new_fasta = cdhit(input.fasta, cdhit_dir; optional_args, cdhit_bin)
        return ASVResult(new_fasta, input.count_table, input.taxonomy)
    end

    function cdhit(project::ProjectCtx, input::ASVResult;
                   optional_args::Union{String,Nothing} = nothing,
                   cdhit_bin = tool_bin("cd_hit_est"))
        config_path   = write_run_config(project)
        cfg           = get(YAML.load_file(config_path), "cdhit", Dict())
        optional_args = isnothing(optional_args) ?
                        get(cfg, "optional_args", "-c 0.9") : optional_args
        cdhit_dir     = joinpath(project.dir, "cdhit")
        new_fasta     = joinpath(cdhit_dir, basename(input.fasta))
        hash_file     = joinpath(cdhit_dir, "config.hash")
        if isfile(new_fasta) &&
           !_section_stale(config_path, "cdhit", hash_file) &&
           mtime(new_fasta) > mtime(input.fasta)
            @info "Skipping cdhit: $new_fasta up to date"
            return ASVResult(new_fasta, input.count_table, input.taxonomy)
        end
        new_fasta = cdhit(input.fasta, cdhit_dir; optional_args, cdhit_bin)
        _write_section_hash(config_path, "cdhit", hash_file)
        pipeline_log(project, "cd-hit-est complete")
        log_written(project, new_fasta)
        return ASVResult(new_fasta, input.count_table, input.taxonomy)
    end
end
