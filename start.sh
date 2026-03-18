#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

## Frontend build

build_frontend() {
  if command -v bun >/dev/null 2>&1; then
    (cd frontend && bun install --frozen-lockfile && bun run build)
    return
  fi

  if [ -x frontend/node_modules/.bin/tsc ] && [ -x frontend/node_modules/.bin/vite ]; then
    (cd frontend && ./node_modules/.bin/tsc && ./node_modules/.bin/vite build)
    return
  fi

  echo "No frontend build toolchain found. Install bun or restore frontend/node_modules." >&2
  exit 1
}

if [ ! -d web/dist ] || [ "${BUILD:-0}" = "1" ]; then
  echo "Building frontend..."
  build_frontend
fi

## Backend

export JULIA_METAMANIFOLD_ROOT="${JULIA_METAMANIFOLD_ROOT:-$(pwd)}"
export JULIA_METAMANIFOLD_PORT="${JULIA_METAMANIFOLD_PORT:-8080}"

JULIA_THREADS="${JULIA_THREADS:-8}"
JULIA_ARGS=(--project=. --threads="$JULIA_THREADS")
if [ -f MetaManifold.so ]; then
  JULIA_ARGS+=(--sysimage MetaManifold.so)
fi

exec julia "${JULIA_ARGS[@]}" src/server/server.jl
