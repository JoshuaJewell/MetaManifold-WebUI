# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

## One-shot stage-hash migration for the r_threads / remote config move
#
# `dada2.taxonomy.multithread` and `dada2.taxonomy.remote` moved to the
# top-level `r_threads` and `remote` blocks. Neither ever changed a single
# taxonomic assignment, but both sat inside the `dada2.taxonomy` section that
# :dada2_assign_taxonomy hashes, so removing them changes that hash and every
# existing run's assign_taxonomy would read as stale - inviting a re-run of the
# most expensive stage in the pipeline to reproduce byte-identical output.
#
# This rewrites each run's stored assign_taxonomy hash to match the config it
# now resolves to, leaving the stage current. It changes no result and no
# checkpoint: only the hash a stage compares itself against.
#
# Run once, after updating config/pipeline.yml:
#
#     julia --project=. -e 'include("scripts/migrate_r_threads.jl"); MigrateRThreads.run_migration("projects")'
#
# Pass `dry_run=true` to list what would be rewritten without touching anything.
module MigrateRThreads

using YAML

using MetaManifold.Config: _write_section_hash, stage_sections

# A run directory is one holding a resolved run_config.yml beside a dada2
# workspace; that pairing is what every stage's hash is computed from.
function _run_dirs(root::String)
    dirs = String[]
    isdir(root) || return dirs
    for (dir, _, files) in walkdir(root)
        "run_config.yml" in files && isdir(joinpath(dir, "dada2", "Checkpoints")) &&
            push!(dirs, dir)
    end
    return sort(dirs)
end

"""
    run_migration(projects_root; dry_run=false) -> NamedTuple

Rewrite the stored assign_taxonomy stage hash of every run under
`projects_root` so the move of `multithread` and `remote` out of
`dada2.taxonomy` does not mark completed work stale.

Only runs that actually completed the stage are touched: one that never
assigned taxonomy has no hash file, and writing one would claim work that was
never done.
"""
function run_migration(projects_root::String; dry_run::Bool=false)
    rewritten = String[]
    skipped   = String[]

    for run_dir in _run_dirs(projects_root)
        config_path = joinpath(run_dir, "run_config.yml")
        hash_file   = joinpath(run_dir, "dada2", "Checkpoints", "assign_taxonomy.hash")
        checkpoint  = joinpath(run_dir, "dada2", "Checkpoints", "checkpoint.RData")

        # No stored hash, or no checkpoint behind it, means the stage never
        # finished here. Leave it stale: that is the truth.
        if !isfile(hash_file) || !isfile(checkpoint)
            push!(skipped, run_dir)
            continue
        end

        push!(rewritten, run_dir)
        dry_run || _write_section_hash(config_path,
                                       stage_sections(:dada2_assign_taxonomy),
                                       hash_file)
    end

    (; rewritten, skipped, dry_run)
end

end # module
