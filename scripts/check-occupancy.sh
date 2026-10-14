#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Gate for the occupancy reading of the vsearch stage (docs/occupancy-2026-10.adoc):
#   1. each Session IR text checks, and occupancy-types' stepper counts its bytes;
#   2. TypedPeak.agda type-checks against the pinned siblings, and each module in
#      docs/occupancy/reject is refused at its own declaration;
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
agda "${flags[@]}" docs/occupancy/TypedPeak.agda >"$logs/TypedPeak.log" 2>&1 ||
  { cat "$logs/TypedPeak.log"; echo "ERROR: TypedPeak.agda does not type-check" >&2; exit 1; }
echo "PASS: TypedPeak.agda"
for file in docs/occupancy/reject/*.agda; do
  name="$(basename -- "$file" .agda)"
  log="$logs/$name.log"
  if agda "${flags[@]}" -i docs/occupancy/reject "$file" >"$log" 2>&1; then
    cat "$log"; echo "ERROR: invalid claim accepted: $name" >&2; failures=$((failures + 1)); continue
  fi
  if ! grep -Eq "/$name\.agda:[0-9]+,[0-9]+" "$log" || ! grep -q 'when checking that' "$log" ||
     grep -Eq 'Not in scope|Failed to find|Unsolved|Parse error' "$log"; then
    cat "$log"; echo "ERROR: $name failed for the wrong reason" >&2; failures=$((failures + 1)); continue
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
  fi
done

if (( failures > 0 )); then
  echo "check-occupancy: $failures failure(s)" >&2
  exit 1
fi
echo "check-occupancy: OK"
