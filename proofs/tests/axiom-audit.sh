#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Audit the proof tree for anything that would let a "proof" pass without
# proving anything.
#
# A type-checker only guarantees what its rules allow.  Every check here is a
# rule we refuse to allow.  Each one exits non-zero on a finding, and each one is
# exercised by `proofs/tests/gate-selftest.sh`, so this script cannot rot into a
# no-op unnoticed.
#
# Checks:
#   1. every module declares `--safe`
#   2. no `postulate` blocks
#   3. no FFI, no `{-# COMPILE`, no `{-# BUILTIN`
#   4. no unsound flags (`--type-in-type`, `--no-positivity-check`,
#      `--no-termination-check`, `--no-universe-polymorphism` off-switches)
#   5. no `trustMe` / `primTrust` / `{-# NO_POSITIVITY_CHECK` pragmas
#   6. no holes (`{!!}`, `?`) left in a checked module
#   7. every `.agda` file is reachable from `All.agda` — a module that is not
#      imported by the gate entry point is not checked at all, which is the
#      easiest way to have a "proved" theorem that nobody ever type-checked

set -euo pipefail

PROOFS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
AGDA_DIR="$PROOFS_DIR/agda"
ENTRY="$AGDA_DIR/MetaManifold/All.agda"

failures=0
fail() { printf 'axiom-audit: FAIL: %s\n' "$*" >&2; failures=$((failures + 1)); }
note() { printf 'axiom-audit: %s\n' "$*"; }

[[ -f "$ENTRY" ]] || { printf 'axiom-audit: FATAL: missing %s\n' "$ENTRY" >&2; exit 1; }

mapfile -t FILES < <(find "$AGDA_DIR/MetaManifold" -name '*.agda' -type f | sort)
[[ ${#FILES[@]} -gt 0 ]] || { printf 'axiom-audit: FATAL: no .agda files found\n' >&2; exit 1; }
note "auditing ${#FILES[@]} module(s)"

for f in "${FILES[@]}"; do
  rel="${f#"$AGDA_DIR"/}"

  # 1. --safe
  if ! grep -q '^{-# OPTIONS .*--safe' "$f"; then
    fail "$rel does not declare --safe"
  fi

  # 2. postulates
  if grep -nE '^[[:space:]]*postulate\b' "$f" >/dev/null; then
    fail "$rel contains a postulate block: $(grep -nE '^[[:space:]]*postulate\b' "$f" | head -3 | tr '\n' ' ')"
  fi

  # 3. FFI
  if grep -nE '\{-#[[:space:]]*(FOREIGN|COMPILE|BUILTIN)' "$f" >/dev/null; then
    fail "$rel uses FOREIGN/COMPILE/BUILTIN"
  fi

  # 4. unsound flags
  if grep -nE '\{-#[[:space:]]*OPTIONS.*(--type-in-type|--no-positivity-check|--no-termination-check|--no-guardedness)' "$f" >/dev/null; then
    fail "$rel enables an unsound flag"
  fi

  # 5. trusted primitives
  if grep -nE '\b(trustMe|primTrust|NO_POSITIVITY_CHECK|NO_TERMINATION_CHECK)\b' "$f" >/dev/null; then
    fail "$rel uses trustMe/primTrust or a check-disabling pragma"
  fi

  # 6. holes
  if grep -nE '\{\![^!]*\!\}|^[^[:space:]].*[[:space:]]\?[[:space:]]*$' "$f" >/dev/null; then
    fail "$rel contains a hole"
  fi
done

# 7. reachability from the gate entry point
#
# Walk the import graph starting at All.agda.  Anything not reached is a module
# the gate never type-checks.
declare -A SEEN=()
queue=("MetaManifold.All")
while [[ ${#queue[@]} -gt 0 ]]; do
  mod="${queue[0]}"
  queue=("${queue[@]:1}")
  [[ -n "${SEEN[$mod]:-}" ]] && continue
  SEEN[$mod]=1
  path="$AGDA_DIR/$(echo "$mod" | tr '.' '/').agda"
  [[ -f "$path" ]] || { fail "imported module $mod has no file at $path"; continue; }
  while read -r dep; do
    queue+=("$dep")
  done < <(grep -oE '^[[:space:]]*(open[[:space:]]+)?import[[:space:]]+[A-Za-z0-9_.]+' "$path" \
             | awk '{print $NF}' | grep '^MetaManifold\.')
done

for f in "${FILES[@]}"; do
  rel="${f#"$AGDA_DIR"/}"
  mod="$(echo "${rel%.agda}" | tr '/' '.')"
  if [[ -z "${SEEN[$mod]:-}" ]]; then
    fail "$rel is not reachable from MetaManifold/All.agda — it is never type-checked"
  fi
done
note "reachability: ${#SEEN[@]} module(s) reachable from MetaManifold.All"

if [[ $failures -gt 0 ]]; then
  printf 'axiom-audit: %d finding(s)\n' "$failures" >&2
  exit 1
fi
note "clean"
