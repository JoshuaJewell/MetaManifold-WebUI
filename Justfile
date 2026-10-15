# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# MetaManifold's task entry (RSR scheme: a thin Justfile at the root; the
# heavy work lives in the commands it names). Every recipe delegates to the
# repository's own documented entry points -- install.sh, start.sh, the test
# runner -- so there is exactly one place where each of those is defined and
# this file never re-implements them.

# The default recipe lists what exists; `just --list` does the same.
default:
    @just --list

# The documented installer: Julia (pinned), R, the external tools, the
# frontend build. Extra args pass through to install.sh (--update,
# --modify, --sysimage).
install *ARGS:
    bash install.sh {{ARGS}}

# The documented server entry point. Serves the UI on port 8080.
serve:
    bash start.sh

# The full unit suite. Name files to run a subset:
#   just test test_install_pins.jl test_server_security.jl
test *FILES:
    julia --project=. test/runtests.jl {{FILES}}

# The pin suite alone: toolchain pins, lockfile agreement, workflow shape,
# and the root-shape allowlist.
pins:
    julia --project=. test/runtests.jl test_install_pins.jl

# The Agda proof tree's audit (no postulates, no holes, everything reached
# from All.agda). The type-check itself runs in CI's Agda proofs job.
proofs:
    bash proofs/tests/axiom-audit.sh

# The frontend's own checks.
frontend-test:
    cd frontend && bun run test

frontend-typecheck:
    cd frontend && bun run typecheck

frontend-build:
    cd frontend && bun run build

# The git-signing agent's state, as the launcher sees it.
agent-status:
    ./launcher.sh --agent status
