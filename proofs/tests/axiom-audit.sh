#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Audit the Agda proof tree for anything that would let a "proof" pass
# without proving anything.  The type-checker only guarantees what its rules
# allow; every check here is a rule this project refuses to relax.
#
# Usage: proofs/tests/axiom-audit.sh [AGDA_ROOT]
#   AGDA_ROOT defaults to proofs/agda next to this script.  Passing another
#   root is how the positive control runs the audit against a scratch copy.
#
# Fails (exit 1) on:
#   1. a module without `--safe` and `--without-K` in its OPTIONS pragma
#   2. a `postulate`, anywhere in code (not only at the start of a line)
#   3. FOREIGN / COMPILE / BUILTIN pragmas
#   4. unsound flags (--type-in-type, --no-positivity-check,
#      --no-termination-check, --allow-unsolved-metas, ...)
#   5. TERMINATING / NON_TERMINATING / NO_POSITIVITY_CHECK /
#      NO_UNIVERSE_CHECK pragmas, trustMe, primTrust
#   6. holes: `{! !}` or a bare `?` token, anywhere in code
#   7. a module under the root that MetaManifold/All.agda does not reach
#      by imports (an unreached module is never type-checked at all)
# Line comments (`-- ...`) are stripped before checks 2 and 6, so prose that
# mentions a postulate is not a finding.  Block comments are not stripped:
# the proof modules do not use them, and a false positive is the safe side.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AGDA_DIR="$(cd "${1:-$SCRIPT_DIR/../agda}" && pwd)"
ENTRY="$AGDA_DIR/MetaManifold/All.agda"
failures=0

# Record one finding on stderr and count it; the audit fails at the end.
fail() {
  printf 'axiom-audit: FAIL: %s\n' "$*" >&2
  failures=$((failures + 1))
}

# Print one progress line on stdout.
note() {
  printf 'axiom-audit: %s\n' "$*"
}

# Abort at once on a condition that makes the audit meaningless.
fatal() {
  printf 'axiom-audit: FATAL: %s\n' "$*" >&2
  exit 1
}

# Print a file with every `-- ...` line comment removed, keeping line numbers.
# A `--` glued to a word (an option such as --safe) is not a comment.
code_only() {
  sed -E 's/(^|[[:space:]])--([[:space:]].*)?$/\1/' "$1"
}

# Run checks 1-6 on one module file, recording findings via `fail`.
audit_file() {
  local f="$1" rel="${1#"$AGDA_DIR"/}" code hits
  code="$(code_only "$f")"

  if ! grep -qE '^\{-#[[:space:]]*OPTIONS.*--safe' "$f"; then
    fail "$rel does not declare --safe"
  fi
  if ! grep -qE '^\{-#[[:space:]]*OPTIONS.*--without-K' "$f"; then
    fail "$rel does not declare --without-K"
  fi

  if hits="$(grep -nE '(^|[^[:alnum:]_])postulate([^[:alnum:]_]|$)' <<<"$code")"; then
    fail "$rel contains a postulate: $(head -3 <<<"$hits" | tr '\n' ' ')"
  fi

  if grep -qE '\{-#[[:space:]]*(FOREIGN|COMPILE|BUILTIN)' "$f"; then
    fail "$rel uses a FOREIGN/COMPILE/BUILTIN pragma"
  fi

  if grep -qE -- '--(type-in-type|no-positivity-check|no-termination-check|allow-unsolved-metas|allow-incomplete-matches|no-guardedness|injective-type-constructors|omega-in-omega|cumulativity)' "$f"; then
    fail "$rel enables an unsound or check-disabling flag"
  fi

  if hits="$(grep -nE '\{-#[[:space:]]*(TERMINATING|NON_TERMINATING|NO_POSITIVITY_CHECK|NO_UNIVERSE_CHECK|NON_COVERING|INJECTIVE)\b|\b(trustMe|primTrust)\b' "$f")"; then
    fail "$rel uses a check-disabling pragma or trusted primitive: $(head -3 <<<"$hits" | tr '\n' ' ')"
  fi

  if hits="$(grep -nE '\{![^!]*!\}|(^|[[:space:](])\?([[:space:])]|$)' <<<"$code")"; then
    fail "$rel contains a hole: $(head -3 <<<"$hits" | tr '\n' ' ')"
  fi
}

# Print, one per line, every MetaManifold.* module reachable from
# MetaManifold.All by imports.  An import with no file is printed as
# `MISSING <module>` (this runs in a subshell, so it cannot call `fail`).
reachable_modules() {
  local -A seen=()
  local -a queue=("MetaManifold.All")
  local mod path dep
  while [[ ${#queue[@]} -gt 0 ]]; do
    mod="${queue[0]}"
    queue=("${queue[@]:1}")
    [[ -n "${seen[$mod]:-}" ]] && continue
    seen[$mod]=1
    printf '%s\n' "$mod"
    path="$AGDA_DIR/${mod//.//}.agda"
    [[ -f "$path" ]] || { printf "MISSING %s\n" "$mod"; continue; }
    while read -r dep; do
      [[ -n "$dep" ]] && queue+=("$dep")
    done < <(code_only "$path" \
               | grep -oE '^[[:space:]]*(open[[:space:]]+)?import[[:space:]]+MetaManifold\.[A-Za-z0-9_.]+' \
               | awk '{print $NF}' || true)
  done
}

[[ -f "$ENTRY" ]] || fatal "missing $ENTRY"
mapfile -t FILES < <(find "$AGDA_DIR" -name '*.agda' -type f | sort)
[[ ${#FILES[@]} -gt 0 ]] || fatal "no .agda files under $AGDA_DIR"
note "auditing ${#FILES[@]} module(s) under $AGDA_DIR"

for f in "${FILES[@]}"; do
  audit_file "$f"
done

REACHED="$(reachable_modules)"
while read -r missing; do
  fail "imported module ${missing#MISSING } has no file"
done < <(grep "^MISSING " <<<"$REACHED" || true)
for f in "${FILES[@]}"; do
  rel="${f#"$AGDA_DIR"/}"
  mod="${rel%.agda}"
  mod="${mod//\//.}"
  if ! grep -qxF "$mod" <<<"$REACHED"; then
    fail "$rel is not reachable from MetaManifold/All.agda, so it is never type-checked"
  fi
done
note "reachability: $(wc -l <<<"$REACHED") module(s) reachable from MetaManifold.All"

if [[ $failures -gt 0 ]]; then
  printf 'axiom-audit: %d finding(s)\n' "$failures" >&2
  exit 1
fi
note "clean"
