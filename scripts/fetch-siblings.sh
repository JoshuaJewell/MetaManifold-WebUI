#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Fetch the proof siblings the occupancy reading (docs/occupancy-2026-10.adoc)
# consumes, each at its pinned commit, into _build/dependencies. Their
# definitions are used, never copied. Re-pin one sibling per commit.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."

# occupancy-types: the HWM Timeline (Agda) and the SessionIR stepper (Julia).
OCCUPANCY_TYPES_COMMIT="${OCCUPANCY_TYPES_COMMIT:-af957280075699c659de0abd9a3b0defade77309}"
# resource-algebra: the HWM grade algebra occupancy-types builds on.
RESOURCE_ALGEBRA_COMMIT="${RESOURCE_ALGEBRA_COMMIT:-cf9235c532712a783a404c359a1055f40134794a}"

fetch() {
  local name="$1" commit="$2" dir="_build/dependencies/$1"
  if [ ! -d "$dir/.git" ]; then
    git init --quiet "$dir"
    git -C "$dir" remote add origin "https://github.com/hyperpolymath/$name.git"
  fi
  git -C "$dir" fetch --quiet --depth 1 origin "$commit"
  git -C "$dir" checkout --quiet --detach "$commit"
  echo "$name at $(git -C "$dir" rev-parse --short HEAD)"
}

mkdir -p _build/dependencies
fetch occupancy-types "$OCCUPANCY_TYPES_COMMIT"
fetch resource-algebra "$RESOURCE_ALGEBRA_COMMIT"
