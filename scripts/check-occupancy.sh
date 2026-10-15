#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Gate for the occupancy reading of the vsearch stage (docs/occupancy-2026-10.adoc):
#   1. each Session IR text checks, and occupancy-types' stepper counts its bytes;
#   2. TypedPeak.agda and GradedStages.agda type-check against the pinned
#      siblings, and each module in docs/occupancy/reject is refused at its own
#      declaration;
#   3. the counted bytes equal the typed peak, order by order.
# Needs julia and agda (2.6.4.x) on PATH. Exit 0 only if all three hold.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
export LC_ALL="${LC_ALL:-C.UTF-8}"

bash scripts/fetch-siblings.sh
deps=_build/dependencies
sir="$deps/occupancy-types/src/julia"
agda --version | head -1
flags=(--safe --without-K --no-libraries --ignore-interfaces --double-check
       -i docs/occupancy -i "$deps/occupancy-types/src" -i "$deps/resource-algebra/src")
logs=_build/check-logs
mkdir -p "$logs"
failures=0

# GitHub renders an `::error` line as a check annotation, and annotations -- unlike
# job logs -- can be read back through the API. A failure that only happens in CI
# must be readable from a machine that cannot download the log, so every failure
# below says what it was and repeats the prover's own last lines: % escapes to
# %25 and newlines to %0A, as the workflow-command syntax requires.
annotate() {  # title file log
  [ -n "${GITHUB_ACTIONS:-}" ] || return 0
  printf '::error file=%s,title=%s::%s\n' "$2" "$1" \
    "$(awk 'NR > 1 { printf "%%0A" } { gsub(/%/, "%25"); printf "%s", $0 }' "$3" |
       tail -c 6000)"
}

# One Agda module that must type-check, with its own log.
check_agda() {  # file label
  if ! agda "${flags[@]}" "$1" >"$logs/$2.log" 2>&1; then
    cat "$logs/$2.log"
    annotate "Agda refused $2" "$1" "$logs/$2.log"
    echo "ERROR: $1 does not type-check" >&2
    exit 1
  fi
  echo "PASS: $2"
}

# Order -> Session IR file and TypedPeak lemma.
orders=(asv-then-otu otu-then-asv concurrent)

echo "== Session IR (occupancy-types v0 stepper)"
declare -A counted
for o in "${orders[@]}"; do
  file="docs/occupancy/vsearch-$o.sir"
  julia --startup-file=no --project="$sir" "$sir/bin/sessionir.jl" check "$file"
  out="$(julia --startup-file=no --project="$sir" "$sir/bin/sessionir.jl" run "$file" --trace)"
  echo "$out"
  counted[$o]="$(sed -n 's/.*peak_live=\([0-9]*\).*/\1/p' <<<"$out" | head -1)"
done

echo "== Typed peak (Agda, HWM Timeline)"
check_agda docs/occupancy/TypedPeak.agda TypedPeak.agda

# The graded reading (fork issue #105): occupancy-types' Reclaim theorems at the
# measured grades, from the same pinned siblings. The recorded run's grades
# against its measurement are checked by scripts/check-occupancy-grades.jl, which
# needs the project's Julia packages and so runs in the disagreement job.
check_agda docs/occupancy/GradedStages.agda GradedStages.agda
for file in docs/occupancy/reject/*.agda; do
  name="$(basename -- "$file" .agda)"
  log="$logs/$name.log"
  if agda "${flags[@]}" -i docs/occupancy/reject "$file" >"$log" 2>&1; then
    cat "$log"
    annotate "Invalid claim accepted: $name" "$file" "$log"
    echo "ERROR: invalid claim accepted: $name" >&2; failures=$((failures + 1)); continue
  fi
  if ! grep -Eq "/$name\.agda:[0-9]+,[0-9]+" "$log" || ! grep -q 'when checking that' "$log" ||
     grep -Eq 'Not in scope|Failed to find|Unsolved|Parse error' "$log"; then
    cat "$log"
    annotate "$name failed for the wrong reason" "$file" "$log"
    echo "ERROR: $name failed for the wrong reason" >&2; failures=$((failures + 1)); continue
  fi
  echo "PASS: rejected $name ($(grep -m1 -E '[0-9]+ != [0-9]+' "$log"))"
done

echo "== Counted bytes against typed peak (KiB)"
printf '%-14s %10s %10s\n' order counted typed
for o in "${orders[@]}"; do
  typed="$(grep -A1 -E "^$o :" docs/occupancy/TypedPeak.agda | sed -n 's/.*≡ \([0-9]*\)$/\1/p' | head -1)"
  printf '%-14s %10s %10s\n' "$o" "${counted[$o]}" "$typed"
  if [ -z "${counted[$o]}" ] || [ -z "$typed" ] || [ "${counted[$o]}" != "$typed" ]; then
    echo "ERROR: $o: counted ${counted[$o]:-?} != typed ${typed:-?}" >&2; failures=$((failures + 1))
    echo "counted ${counted[$o]:-?} KiB != typed ${typed:-?} KiB" >"$logs/counted-$o.log"
    annotate "Counted bytes != typed peak ($o)" docs/occupancy/TypedPeak.agda "$logs/counted-$o.log"
  fi
done

if (( failures > 0 )); then
  echo "check-occupancy: $failures failure(s)" >&2
  exit 1
fi
echo "check-occupancy: OK"
