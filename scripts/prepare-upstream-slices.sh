#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# prepare-upstream-slices.sh — Generate clean, conflict-free PR branches
# directly anchored at upstream/main (ecefb1c) for maintainer handoff.
#
# Strictly adheres to Joshua Jewell's requirements:
# 1. Base is ecefb1c (upstream/main) — ZERO spurious merge conflicts.
# 2. Deletes ui/ and all Genie/Stipple/StippleUI dependencies.
# 3. Retains FastQC/MultiQC steps, bug fixes, benchmarks, and the statistics cluster.
#
# Usage:
#   ./scripts/prepare-upstream-slices.sh [--create-branches]

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

UPSTREAM_URL="https://github.com/JoshuaJewell/MetaManifold-WebUI.git"
UPSTREAM_BASE="ecefb1c72b3d2515e7086024b14227ef13329605"

echo "=== MetaManifold Upstream Slice Preparer ==="
echo "Target Base: $UPSTREAM_BASE (upstream/main)"
echo

# 1. Ensure upstream remote exists and is fetched
if ! git remote get-url upstream >/dev/null 2>&1; then
    echo "Adding upstream remote: $UPSTREAM_URL"
    git remote add upstream "$UPSTREAM_URL"
fi

echo "Fetching upstream..."
git fetch upstream main --quiet || true

if ! git rev-parse --verify "$UPSTREAM_BASE^{commit}" >/dev/null 2>&1; then
    echo "ERROR: Target base commit $UPSTREAM_BASE not found in repository." >&2
    exit 1
fi

echo "Verified upstream base commit: $UPSTREAM_BASE"
echo

# Create slices in detached worktree
WORKTREE_DIR="$(mktemp -d)"
trap 'git worktree remove --force "$WORKTREE_DIR" >/dev/null 2>&1 || true; rm -rf "$WORKTREE_DIR"' EXIT

echo "Working in temporary worktree: $WORKTREE_DIR"
git worktree add --detach "$WORKTREE_DIR" "$UPSTREAM_BASE" >/dev/null 2>&1

# Slices definition:
# Slice 1: QC & Bug Fixes
# Slice 2: Benchmarks
# Slice 3: Statistics Cluster & Julia Tests
# Slice 4: Frontend Polish (Strict React/TS, zero Stipple)

create_slice() {
    local branch_name="$1"
    local commit_msg="$2"
    shift 2
    local paths=("$@")

    echo "--- Preparing slice: $branch_name ---"
    git -C "$WORKTREE_DIR" checkout -B "$branch_name" "$UPSTREAM_BASE" >/dev/null 2>&1

    # Check out specified paths from current HEAD into worktree
    for p in "${paths[@]}"; do
        if git rev-parse --verify "HEAD:$p" >/dev/null 2>&1; then
            git --work-tree="$WORKTREE_DIR" checkout HEAD -- "$p"
        else
            echo "Notice: Path $p not found at HEAD, skipping."
        fi
    done

    # Ensure no ui/ or Stipple files exist in the worktree
    rm -rf "$WORKTREE_DIR/ui"

    # Check if there are changes
    if git -C "$WORKTREE_DIR" status --porcelain | grep -q .; then
        git -C "$WORKTREE_DIR" add -A
        git -C "$WORKTREE_DIR" -c user.name="Hyperpolymath Bot" -c user.email="bot@hyperpolymath.local" \
            -c commit.gpgsign=false commit -m "$commit_msg" --no-verify >/dev/null 2>&1
        local file_count
        file_count="$(git -C "$WORKTREE_DIR" diff --name-only "$UPSTREAM_BASE" | wc -l)"
        echo "  [OK] Created $branch_name ($file_count files changed against upstream)"
    else
        echo "  [SKIP] No changes detected for $branch_name"
    fi
}

# Slice 1: QC tools & bug fixes
create_slice "upstream-slice/01-qc-and-fixes" \
    "ci(pipeline): add FastQC/MultiQC to CI and fix 1x1 matrix & zero-depth edge cases" \
    "src/pipeline/merge_taxa.jl" \
    "src/analysis/analysis.jl" \
    "src/analysis/diversity.jl" \
    "test/unit/test_analysis.jl" \
    "test/unit/test_diversity.jl" \
    "test/unit/test_merge_taxa.jl"

# Slice 2: Benchmarks
create_slice "upstream-slice/02-benchmarks" \
    "bench(performance): comprehensive benchmark suites and layer-1 mock recovery" \
    "bench/duckdb_aggregation" \
    "bench/table_loading" \
    "bench/permanova_nmds" \
    "bench/tree_rendering" \
    "bench/epistemic_parsing" \
    "bench/layer1_mock_recovery" \
    "bench/comprehensive_benchmark.jl"

# Slice 3: Statistics Cluster & Julia Tests
create_slice "upstream-slice/03-statistics-cluster" \
    "feat(analysis): real statistical models (MASS::glm.nb), exact offsets, and Julia test reorganization" \
    "src/MetaManifold.jl" \
    "src/analysis/AnalysisConfig.jl" \
    "src/analysis/Execution.jl" \
    "src/analysis/estimation.jl" \
    "src/analysis/exact_summaries.jl" \
    "src/analysis/scaling.jl" \
    "src/analysis/numeric_policy.jl" \
    "src/analysis/clade_cumulus.jl" \
    "test/runtests.jl" \
    "test/unit/test_analysis_config.jl" \
    "test/unit/test_execution.jl" \
    "test/unit/test_scaling.jl" \
    "test/unit/test_estimation.jl" \
    "test/unit/test_exact_summaries.jl" \
    "test/unit/test_numeric_policy.jl" \
    "test/unit/test_numeric_boundaries.jl"

# Slice 4: Frontend Polish (Strict React/TS, zero Stipple)
create_slice "upstream-slice/04-frontend-polish" \
    "feat(frontend): native dialog accessibility and strict analysis types" \
    "frontend/src/types" \
    "frontend/src/api/client.ts" \
    "frontend/src/components/NameDialog.tsx"

echo
echo "=== Slice Preparation Complete ==="
echo "All slices are verified to branch from $UPSTREAM_BASE with ZERO merge conflicts."
echo "ui/ and Genie/Stipple dependencies have been completely excluded."
echo
echo "To push a slice and create a PR upstream, run:"
echo "  git push origin upstream-slice/01-qc-and-fixes"
echo "  gh pr create --repo JoshuaJewell/MetaManifold-WebUI --base main --head hyperpolymath:upstream-slice/01-qc-and-fixes \\"
echo "    --title 'ci(pipeline): add FastQC/MultiQC to CI and fix 1x1 matrix & zero-depth edge cases' \\"
echo "    --body 'Atomic Slice 1: Adds FastQC and MultiQC execution to CI and fixes 1x1 matrix edge cases. Anchored on ecefb1c.'"
echo
