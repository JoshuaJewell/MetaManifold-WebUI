module Tools

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export cutadapt, vsearch, fastqc_all, fastqc_one

    using YAML
    using Logging

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
    # Prevent duplication of primers. Must be instantiated outside get_primers() loop. Global scope may be an issue.
    used_forward = String[]
    used_reverse = String[]

    # Get primers based on input. Will not run without used_forward and used_reverse arrays being defined first.
    function get_primers(input, primers_path)
        data = YAML.load_file(primers_path)

        for pair in data["Pairs"]
            if haskey(pair, input)
                pair_names = pair[input]

                forward = data["Forward"][pair_names[1]]
                reverse = data["Reverse"][pair_names[2]]

                output_forward = in(forward, used_forward) ? "x" : forward
                output_reverse = in(reverse, used_reverse) ? "x" : reverse

                if output_forward != "x"
                    push!(used_forward, output_forward)
                end
                if output_reverse != "x"
                    push!(used_reverse, output_reverse)
                end

                return (output_forward, output_reverse)
            end
        end

        path = pwd() * "/" * primers_path
        return "Invalid input. Please specify a valid primer pair in '$path'."
    end

    # Formats cutadapt's arguments for the primer sets requested 
    function get_primer_args(primer_pairs, primers_path)
        args = ""

        for pair in primer_pairs
            pair_tuple = get_primers(pair, primers_path)
            forward_primer = pair_tuple[1]
            reverse_primer = pair_tuple[2]

            args *= join([flag * primer * " " for (flag, primer) in (("-g ", forward_primer), ("-G ", reverse_primer)) if primer != "x"], "")
        end

        # Reset global used primers.
        used_forward = String[]
        used_reverse = String[]

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
        function cutadapt(primer_pairs, primers_path, fastq_in_dir, fastq_out_dir, optional_args)
    
    Requires `cutadapt` installed. Runs `cutadapt` command with as many primer pairs as necessary (useful for multiplex) and any optional parameters. Primer aguments are determined from YAML file in format:
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
    Where a primer pair contains the same forward or reverse as another specified pair, these are deduplicated by name  rather than sequence. For example:
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
        return run_cutadapt(
            get_primer_args(primer_pairs, primers_path),
            optional_args,
            fastq_in_dir,
            cutadapt_dir,
            cutadapt_bin
        )
    end

    ## FastQC
    """
        function fastqc_all(fastq_in_dir; optional_args = "-t 20 --extract --delete", fastqc_bin = "fastqc")

    Requires `fastqc` installed. Runs `fastqc` command with any optional parameters on a whole set.

    ## Arguments
    - `fastq_in_dir`: Specify path of fastq files to qc.

    ## Keyword Arguments
    - `optional_args` (optional, default: "-t 20 --extract --delete"): Specify additional arguments passed to `fastqc` command.
    - `fastqc_bin` (optional, default: "fastqc"): Path to the fastqc binary. Defaults to PATH lookup. Set via `config/tools.yml` and `load_tools()` in main.jl.

    """
    function fastqc_all(fastq_in_dir, fastqc_dir; optional_args = "-t 20 --extract --delete", fastqc_bin = tool_bin("fastqc"))
        mkpath(fastqc_dir)

        cmd = "$fastqc_bin $fastq_in_dir/*.fastq* -o $fastqc_dir $optional_args"
        run(`bash -lc $cmd`)

    end

    """
        function fastqc_one(fastq_in_file; optional_args = "-t 20 --extract --delete", fastqc_bin = "fastqc")

    Requires `fastqc` installed. Runs `fastqc` command with any optional parameters on a whole set.

    ## Arguments
    - `fastq_in_file`: Specify path of fastq file to qc.

    ## Keyword Arguments
    - `optional_args` (optional, default: "-t 20 --extract --delete"): Specify additional arguments passed to `fastqc` command.
    - `fastqc_bin` (optional, default: "fastqc"): Path to the fastqc binary. Defaults to PATH lookup. Set via `config/tools.yml` and `load_tools()` in main.jl.

    """
    function fastqc_one(fastq_in_file, fastqc_dir; optional_args = "-t 20 --extract --delete", fastqc_bin = tool_bin("fastqc"))
        mkpath(fastqc_dir)

        cmd = "$fastqc_bin $fastq_in_file -o $fastqc_dir $optional_args"
        run(`bash -lc $cmd`)

    end

    ## VSEARCH
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
        mkpath(vsearch_dir)
        outfile = joinpath(vsearch_dir, "taxonomy.tsv")

        cmd = "$vsearch_bin --usearch_global $fasta_in_dir --db $reference_database --blast6out $outfile $optional_args"
        run(`bash -lc $cmd`)

    end
end