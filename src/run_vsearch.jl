module VSEARCH

export vsearch

    using Logging

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
    function vsearch(fasta_in_dir, reference_database, vsearch_dir; optional_args = "--id 0.75 --query_cov 0.8", vsearch_bin = "vsearch")
        mkpath(vsearch_dir)
        outfile = joinpath(vsearch_dir, "taxonomy.tsv")

        cmd = "$vsearch_bin --usearch_global $fasta_in_dir --db $reference_database --blast6out $outfile $optional_args"
        run(`bash -lc $cmd`)

    end

end