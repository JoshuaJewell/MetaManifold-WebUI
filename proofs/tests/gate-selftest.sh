#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Prove that the proof gate can fail.
#
# A gate that has never been observed to reject anything is not evidence.  This
# script takes a throwaway copy of the proof tree, breaks it in one specific way
# at a time, and requires the gate to reject it.  If any mutation is *accepted*,
# or if a mutation failed to apply at all, this script exits non-zero.
#
# It runs in CI.  A gate whose self-test is skipped is a gate nobody can trust,
# so there is no flag to turn this off.
#
# Usage: proofs/tests/gate-selftest.sh

set -uo pipefail

PROOFS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Resolve Agda the same way bootstrap.sh does, so `just proofs-selftest` works
# after `just proofs-bootstrap` with nothing else on PATH.
VENDOR="${PROOFS_VENDOR:-$PROOFS_DIR/.vendor}"
if [[ -n "${AGDA_BIN:-}" ]]; then
  AGDA="$AGDA_BIN"
elif command -v agda >/dev/null 2>&1; then
  AGDA="$(command -v agda)"
elif [[ -x "$VENDOR/venv/bin/agda" ]]; then
  AGDA="$VENDOR/venv/bin/agda"
else
  printf 'gate-selftest: FATAL: agda not found (run proofs/bootstrap.sh --bootstrap, or set AGDA_BIN)\n' >&2
  exit 1
fi
[[ -x "$AGDA" ]] || { printf 'gate-selftest: FATAL: %s is not executable\n' "$AGDA" >&2; exit 1; }

STDLIB_LIB=""
for cand in "$VENDOR/agda-stdlib/standard-library.agda-lib"; do
  [[ -f "$cand" ]] && { STDLIB_LIB="$cand"; break; }
done
[[ -n "$STDLIB_LIB" ]] || { printf 'gate-selftest: FATAL: standard library not found\n' >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cp -r "$PROOFS_DIR/agda" "$WORK/agda"
cp -r "$PROOFS_DIR/tests" "$WORK/tests"
cp "$PROOFS_DIR/agda/metamanifold.agda-lib" "$WORK/agda/" 2>/dev/null || true

# Point Agda at the *copy*: without this the library file would resolve
# `MetaManifold.All` back to the pristine tree and every mutation below would be
# checked against unmutated sources — a self-test that tests nothing.
printf '%s\n%s\n' "$WORK/agda/metamanifold.agda-lib" "$STDLIB_LIB" > "$WORK/libraries"

run_gate() {
  ( cd "$WORK/agda" && "$AGDA" --library-file="$WORK/libraries" --safe --without-K MetaManifold/All.agda ) \
    >"$WORK/out.txt" 2>&1
}
run_audit() { "$WORK/tests/axiom-audit.sh" >"$WORK/out.txt" 2>&1; }

total=0
failed=0

# expect_reject <name> <runner> <file> <mutation-script>
expect_reject() {
  local name="$1" runner="$2" target="$3" mutator="$4"
  local file="$WORK/agda/$target"
  total=$((total + 1))

  # Restore a pristine copy of the target before mutating.
  cp "$PROOFS_DIR/agda/$target" "$file"
  local before after
  before="$(sha256sum "$file" | awk '{print $1}')"
  bash -c "$mutator" mutator "$file"
  after="$(sha256sum "$file" | awk '{print $1}')"
  if [[ "$before" == "$after" ]]; then
    printf 'gate-selftest: FAIL  %-42s mutation did not apply (stale pattern?)\n' "$name"
    failed=$((failed + 1))
    return
  fi

  if $runner; then
    printf 'gate-selftest: FAIL  %-42s gate ACCEPTED a broken proof\n' "$name"
    sed 's/^/                     | /' "$WORK/out.txt" | head -6
    failed=$((failed + 1))
  else
    printf 'gate-selftest: ok    %-42s rejected\n' "$name"
  fi
  cp "$PROOFS_DIR/agda/$target" "$file"
}

# --- the type-checker must reject wrong mathematics -------------------------
#
# Each mutator is a shell snippet receiving the target file as $1.  If a snippet
# stops matching (because the proof was rewritten), `expect_reject` reports the
# mutation as unapplied rather than silently passing.

expect_reject "rounding: wrong known-answer digit" run_gate "MetaManifold/DecimalRounding.agda" \
  'sed -i "s|+ 402 ℤ.<? + 403|+ 402 ℤ.<? + 402|" "$1"'

expect_reject "rounding: tie forced the wrong way" run_gate "MetaManifold/DecimalRounding.agda" \
  'sed -i "s|half-ties-round-up : IsRoundHalfUp (+ 1) 1 0 (+ 1)|half-ties-round-up : IsRoundHalfUp (+ 1) 1 0 (+ 0)|" "$1"'

expect_reject "bh: monotonicity reversed" run_gate "MetaManifold/BenjaminiHochberg.agda" \
  'sed -i "s|  bhScale M j n d ≤ℚ bhScale M′ j n d|  bhScale M′ j n d ≤ℚ bhScale M j n d|" "$1"'

expect_reject "permutation: plus-one removed" run_gate "MetaManifold/PermutationTest.agda" \
  'sed -i "s|pValue b B = mkℚᵘ (+\[1+ b \]) B|pValue b B = mkℚᵘ (+ b) B|" "$1"'

expect_reject "proportions: zero total silently zero" run_gate "MetaManifold/Proportions.agda" \
  'sed -i "s#^\.\.\. | yes _ = refused zeroTotal#... | yes _ = value 0ℚᵘ#" "$1"'

expect_reject "exact counts: overflow no longer refused" run_gate "MetaManifold/ExactCounts.agda" \
  'sed -i "s|checkedAdd-refuses-exactly-when-it-must :|checkedAdd-refuses-exactly-when-it-must-DISABLED :|" "$1"'

# --- the axiom audit must reject an unchecked or unsound module -------------

expect_reject "audit: module dropped from the gate entry" run_audit "MetaManifold/All.agda" \
  'sed -i "/import MetaManifold.DecimalRounding/d" "$1"'

expect_reject "audit: postulate injected" run_audit "MetaManifold/DecimalRounding.agda" \
  'printf "postulate cheat : ∀ {A : Set} → A\n" >> "$1"'

expect_reject "audit: --safe removed" run_audit "MetaManifold/Proportions.agda" \
  'sed -i "s|--without-K --safe|--without-K|" "$1"'

# --- the gate must still accept the pristine tree --------------------------

total=$((total + 1))
cp "$PROOFS_DIR/agda/MetaManifold/All.agda" "$WORK/agda/MetaManifold/All.agda"
if run_gate && run_audit; then
  printf 'gate-selftest: ok    %-42s accepted\n' "pristine tree"
else
  printf 'gate-selftest: FAIL  %-42s gate REJECTED the real proofs\n' "pristine tree"
  sed 's/^/                     | /' "$WORK/out.txt" | head -12
  failed=$((failed + 1))
fi

printf 'gate-selftest: %d/%d controls behaved correctly\n' "$((total - failed))" "$total"
[[ $failed -eq 0 ]] || exit 1
