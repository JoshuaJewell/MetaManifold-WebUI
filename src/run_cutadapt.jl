module Cutadapt

export cutadapt

    using YAML
    using TimeZones
    using Logging

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

        # Messy filesystem stuff
        time = chop("$(now(localzone()))", tail = 13)
        fastq_out_dir = cutadapt_dir * "$time/"
        log_dir = fastq_out_dir * "/logs/"
        stats_path = joinpath(log_dir, "cutadapt_primer_trimming_stats.txt")
        summary_path = joinpath(log_dir, "cutadapt_trimmed_percentage.txt")
        stats_basename = basename(stats_path)
        summary_basename = basename(summary_path)

        isdir(log_dir) || mkpath(log_dir)

        # Filter .fastq.gz files
        for f in filter(x->occursin(r"_R1.*fastq\.gz$", x), readdir(fastq_in_dir))
            push!(samples, split(f, '_')[1])
        end

        @info("cutadapt running at $time with arguments: $primer_args $optional_args.")

        nsamples = length(samples)
        for (i, sample) in enumerate(samples)
            inputR1 = fastq_in_dir * sample * "_*_L001_R1_001.fastq.gz"
            inputR2 = fastq_in_dir * sample * "_*_L001_R2_001.fastq.gz"

            outputR1 = fastq_out_dir * sample * "_R1_trimmed.fastq.gz"
            outputR2 = fastq_out_dir * sample * "_R2_trimmed.fastq.gz"

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

        @info("cutadapt complete. Output available in $fastq_out_dir.")
    end

    """
        function cutadapt(primer_pairs, primers_path, fastq_in_dir, fastq_out_dir, optional_args)
    
    Requires `cutadapt` installed and in PATH. Runs `cutadapt` command with as many primer pairs as necessary (useful for multiplex) and any optional parameters. Primer aguments are determined from YAML file in format:
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
        cutadapt_bin  = "cutadapt"
    )
        run_cutadapt(
            get_primer_args(primer_pairs, primers_path),
            optional_args,
            fastq_in_dir,
            cutadapt_dir,
            cutadapt_bin
        )
    end
end