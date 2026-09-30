#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

## Frontend build
# install.sh puts the pinned bun in bin/; a bun on PATH is the fallback.
BUN=""
if [ -x bin/bun ]; then
  BUN="$(pwd)/bin/bun"
elif command -v bun >/dev/null 2>&1; then
  BUN="$(command -v bun)"
fi

frontend_stale() {
  [ -f web/dist/index.html ] || return 0
  [ -n "$(find frontend/src frontend/public frontend/index.html frontend/package.json \
            frontend/bun.lock frontend/vite.config.ts frontend/tsconfig.json \
            -newer web/dist/index.html -print -quit 2>/dev/null)" ]
}

if [ "${BUILD:-0}" = "1" ] || frontend_stale; then
  if [ -z "$BUN" ]; then
    echo "The frontend needs building and bun was not found. Run ./install.sh first." >&2
    exit 1
  fi
  echo "Building frontend..."
  (cd frontend && "$BUN" install --frozen-lockfile && "$BUN" run build)
fi

## Backend

# install.sh unpacks the pinned release into bin/julia when no julia of the
# pinned version was available, so one there takes precedence.
if [ -x bin/julia/bin/julia ]; then
  PATH="$(pwd)/bin/julia/bin:$PATH"
fi

# Fall back to juliaup's shim dir if julia isn't already on PATH (e.g. a fresh
# shell that never sourced ~/.juliaup/env after install.sh ran).
if ! command -v julia >/dev/null 2>&1; then
  if [ -f "$HOME/.juliaup/env" ]; then
    # shellcheck disable=SC1091
    . "$HOME/.juliaup/env"
  fi
  [ -x "$HOME/.juliaup/bin/julia" ] && PATH="$HOME/.juliaup/bin:$PATH"
fi

if ! command -v julia >/dev/null 2>&1; then
  echo "julia not found on PATH." >&2
  echo "Run ./install.sh, or add juliaup to your PATH:" >&2
  echo '  export PATH="$HOME/.juliaup/bin:$PATH"' >&2
  exit 1
fi

export JULIA_METAMANIFOLD_ROOT="${JULIA_METAMANIFOLD_ROOT:-$(pwd)}"
export JULIA_METAMANIFOLD_PORT="${JULIA_METAMANIFOLD_PORT:-8080}"

JULIA_THREADS="${JULIA_THREADS:-8}"
JULIA_ARGS=(--project=. --threads="$JULIA_THREADS")
for sysimage in MetaManifold.so MetaManifold.dylib; do
  if [ -f "$sysimage" ]; then
    JULIA_ARGS+=(--sysimage "$sysimage")
    break
  fi
done

exec julia "${JULIA_ARGS[@]}" scripts/serve.jl
