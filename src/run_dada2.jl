module OldDADA2

# This is now obsolete, I guess...how sad, it was quite cute.

export dada2

    using RCall

    """
        dada2(dada2_config_dir)

    Most adorable little wrapper that simply passes the dada2 config to R to run through its pipeline. Requires dada2.r in src.

    ## Arguments:
    - `dada2_config_dir`: Path to DADA2 YAML config.
    """
    function dada2(dada2_config_dir)
        R"system(paste('Rscript', './src/dada2.r', $dada2_config_dir))"
    end

end