module DuckDBStore

# © 2026 Joshua Benjamin Jewell. All rights reserved.
#
# This module is licensed under the GNU Affero General Public License version 3 (AGPLv3).

export load_results_db, with_results_db, with_results_db_write

using DuckDB, Logging

"""
    load_results_db(merge_dir; swarm_dir=nothing)

Create/recreate a DuckDB database at `merge_dir/results.duckdb` by loading all
CSVs in the merge directory. If `swarm_dir` is provided and contains
`cluster_membership.csv`, that table is loaded too.
"""
function load_results_db(merge_dir::String; swarm_dir::Union{String,Nothing}=nothing)
    db_path = joinpath(merge_dir, "results.duckdb")
    # Remove stale DB and WAL to avoid schema conflicts on re-runs
    isfile(db_path) && rm(db_path)
    isfile(db_path * ".wal") && rm(db_path * ".wal")

    db  = DuckDB.DB(db_path)
    con = DBInterface.connect(db)
    try
        _esc_id(s) = "\"" * replace(s, "\"" => "\"\"") * "\""
        _esc_str(s) = "'" * replace(s, "'" => "''") * "'"

        for f in readdir(merge_dir)
            endswith(f, ".csv") || continue
            table_name = splitext(f)[1]
            csv_path   = joinpath(merge_dir, f)
            DBInterface.execute(con,
                "CREATE TABLE $(_esc_id(table_name)) AS SELECT * FROM read_csv_auto($(_esc_str(csv_path)), auto_detect=true)")
            @info "DuckDB: loaded $f as table '$table_name'"
        end

        if !isnothing(swarm_dir)
            membership_csv = joinpath(swarm_dir, "cluster_membership.csv")
            if isfile(membership_csv)
                DBInterface.execute(con,
                    "CREATE TABLE cluster_membership AS SELECT * FROM read_csv_auto($(_esc_str(membership_csv)), auto_detect=true)")
                @info "DuckDB: loaded cluster_membership from $membership_csv"
            end
        end
    finally
        DBInterface.close!(con)
        close(db)
    end

    @info "DuckDB: database written to $db_path"
    return db_path
end

"""
    with_results_db(f, merge_dir)

Open the DuckDB database read-only and pass the connection to `f`.
Ensures the connection and DB are closed after use.
"""
function with_results_db(f::Function, merge_dir::String)
    db_path = joinpath(merge_dir, "results.duckdb")
    isfile(db_path) || error("DuckDB file not found: $db_path")
    db  = DuckDB.DB(db_path; readonly=true)
    con = DBInterface.connect(db)
    try
        return f(con)
    finally
        DBInterface.close!(con)
        close(db)
    end
end

"""
    with_results_db_write(f, merge_dir)

Open the DuckDB database in read-write mode and pass the connection to `f`.
Use this for mutations (CREATE TABLE, DROP TABLE) only; prefer `with_results_db`
for all read operations to avoid unintentional writes.
"""
function with_results_db_write(f::Function, merge_dir::String)
    db_path = joinpath(merge_dir, "results.duckdb")
    isfile(db_path) || error("DuckDB file not found: $db_path")
    db  = DuckDB.DB(db_path)
    con = DBInterface.connect(db)
    try
        return f(con)
    finally
        DBInterface.close!(con)
        close(db)
    end
end

end
