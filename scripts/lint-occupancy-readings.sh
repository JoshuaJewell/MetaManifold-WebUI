#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Fork issue #105: a task's recorded `peak_rss_kib` is a reading, not an
# occupancy. What a stage occupied is determined only by what the next stage
# inherits from it — the drift — which is why the probe grades a whole window and
# why src/core/memory_probe.jl's helpers (`graded`, `predicted_peak_rss_kib`,
# `cost_bound_rss_kib`, `reclaimed_kib`, `task_grades`) exist. Code outside the
# probe that reaches for the key itself is comparing measured figures without
# that context; on the recorded run that reading is off by 9180 KiB for
# `asv_then_otu` and 9204 KiB for `otu_then_asv` (the reclaimed suffix).
#
# So: production code may name the key only in src/core/memory_probe.jl, where
# the record is built and graded. Tests may read it (they assert the recorded
# numbers directly, which is fine — they are evidence, not presentation).
#
# Run from the repository root:  bash scripts/lint-occupancy-readings.sh
set -u -o pipefail

if [[ ! -f src/core/memory_probe.jl ]]; then
  echo "lint-occupancy-readings: run from the repository root" >&2
  exit 2
fi

failures=0
while IFS= read -r hit; do
  [[ -z "$hit" ]] && continue
  echo "ERROR: ungraded occupancy reading: $hit" >&2
  failures=$((failures + 1))
done < <(grep -rn --include='*.jl' '"peak_rss_kib"' src/ |
         grep -v '^src/core/memory_probe\.jl:' || true)

if (( failures > 0 )); then
  echo "lint-occupancy-readings: $failures ungraded reading(s); use the grading helpers in" >&2
  echo "MemoryProbe (see docs/occupancy-2026-10.adoc, 'Reading the numbers') instead." >&2
  exit 1
fi

echo "lint-occupancy-readings: OK — every peak_rss_kib reading in src/ goes through the probe's grades"
