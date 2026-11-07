# Writes a small DADA2 reference for CI, cut from the verified PR2 download.
#
# assignTaxonomy holds about 262 KB per distinct lineage. PR2 5.1.1 has 54,122 of them,
# about 14.6 GiB, which a hosted runner cannot spare. An evenly spaced 2,000 needs about 0.5 GiB,
# so CI can run the stage instead of skipping it. The same input gives the same subset every time.

root = joinpath(@__DIR__, "..", "..")
include(joinpath(root, "src/core/types.jl"))
include(joinpath(root, "src/core/databases.jl"))
using .Databases

const SUBSET_SIZE = 2000

function read_records(path)
    records = Pair{String,String}[]
    header, seq = "", IOBuffer()
    for line in eachline(`gzip -dc $path`)
        if startswith(line, ">")
            isempty(header) || push!(records, header => String(take!(seq)))
            header = line
        else
            print(seq, line)
        end
    end
    isempty(header) || push!(records, header => String(take!(seq)))
    return records
end

function main()
    dbs     = Databases.ensure_databases(joinpath(root, "config/ci/databases.yml"))
    records = read_records(dbs["pr2_dada2"])
    # First record of each lineage, in file order, then an even spread across the file.
    seen     = Set{String}()
    distinct = [r for r in records if !(first(r) in seen) && (push!(seen, first(r)); true)]
    step     = max(1, length(distinct) ÷ SUBSET_SIZE)
    picked   = first(distinct[1:step:end], SUBSET_SIZE)
    out      = joinpath(root, "databases", "pr2_ci_subset_dada2.fasta")
    open(out, "w") do io
        for (header, seq) in picked
            println(io, header)
            println(io, seq)
        end
    end
    @info "Wrote $(length(picked)) of $(length(distinct)) lineages to $out"
end

main()
