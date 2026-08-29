#!/usr/bin/env bash
# Run the determinism suite in two fresh Julia processes with different thread
# counts and require identical output fingerprints from each.
#
#   test/determinism/cross_process.sh [julia-binary]
set -euo pipefail
JULIA="${1:-julia}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$(mktemp -d)"
cat > "$OUT/run.jl" <<'JL'
using Test, MetaManifold, CSV, DataFrames, JSON3, YAML, DuckDB, DBInterface
using MetaManifold.PipelineTypes, MetaManifold.Config, MetaManifold.Categories, MetaManifold.CompositionLibrary
include(joinpath(ENV["MM_ROOT"], "test", "unit", "test_determinism.jl"))
JL
for t in 1 8; do
  MM_ROOT="$ROOT" MM_DET_DUMP="$OUT/threads$t.tsv" "$JULIA" --project="$ROOT" -t "$t" "$OUT/run.jl"
done
if cmp -s "$OUT/threads1.tsv" "$OUT/threads8.tsv"; then
  echo "cross-process: $(wc -l < "$OUT/threads1.tsv") outputs identical between -t1 and -t8 processes"
else
  echo "cross-process: outputs DIFFER"; diff "$OUT/threads1.tsv" "$OUT/threads8.tsv"; exit 1
fi
