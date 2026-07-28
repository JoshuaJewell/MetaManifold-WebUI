#!/usr/bin/env bash
#
# Bootstrap installer for MetabarcodingPipeline
#
# Checks for Julia and R, installs Julia via juliaup if missing,
# then hands off to install.jl for all further dependency setup.
#
# Usage:
#   bash install.sh [--update] [--modify] [--sysimage]
#
# Options:
#   --update    Refuse if tracked files have uncommitted changes; otherwise abort
#               any in-progress merge, sync tracked files to origin/main,
#               clean generated build artifacts, then re-assert the pinned tool
#               versions in bin/ (see config/defaults/tool_versions.yml). This
#               re-fetches the pins; it does not advance them.
#   --modify    Revisit configured tool paths instead of silently reusing them
#   --sysimage  Pass through to install.jl to build the Julia sysimage

set -euo pipefail

UPDATE_MODE=false
for arg in "$@"; do
    if [ "$arg" = "--update" ]; then
        UPDATE_MODE=true
        break
    fi
done

# OS detection

OS="$(uname -s)"
case "$OS" in
    Linux*)  OS_TYPE="Linux"  ;;
    Darwin*) OS_TYPE="macOS"  ;;
    *)
        echo "Unsupported OS: $OS"
        echo "This installer supports Linux and macOS only."
        exit 1
        ;;
esac

echo "Detected OS: $OS_TYPE"
echo ""

update_checkout() {
    if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        echo "--update requested, but this directory is not a Git worktree."
        echo "Skipping checkout sync."
        echo ""
        return
    fi

    if ! git remote get-url origin >/dev/null 2>&1; then
        echo "--update requested, but no 'origin' remote is configured."
        echo "Skipping checkout sync."
        echo ""
        return
    fi

    if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
        echo "--update would discard uncommitted changes to tracked files:"
        git status --short --untracked-files=no
        echo "Commit or stash them, then re-run: bash install.sh --update"
        exit 1
    fi

    echo "Syncing checkout to origin/main..."
    git merge --abort >/dev/null 2>&1 || true
    git fetch origin
    git reset --hard origin/main
    echo ""
}

if [ "$UPDATE_MODE" = true ]; then
    update_checkout
fi

# Julia check and install

JULIAUP_BIN="$HOME/.juliaup/bin"

# Locate julia, pulling juliaup's shim dir onto PATH if that is all that's
# missing. Returns 0 with julia callable, or 1 if it genuinely can't be found.
locate_julia() {
    command -v julia &>/dev/null && return 0

    if [ -f "$HOME/.juliaup/env" ]; then
        # shellcheck disable=SC1091
        source "$HOME/.juliaup/env"
    fi

    if [ -x "$JULIAUP_BIN/julia" ]; then
        # Standard juliaup layout: installed, but its bin dir isn't on PATH.
        # This is the usual reason install.sh appears to work but a later
        # start.sh in a fresh shell cannot find julia.
        export PATH="$JULIAUP_BIN:$PATH"
    elif command -v juliaup &>/dev/null; then
        # juliaup from a distro package / brew / snap: different shim location.
        juliaup add release >/dev/null 2>&1 || true
        export PATH="$(CDPATH= cd -- "$(dirname -- "$(command -v juliaup)")" && pwd):$PATH"
    fi

    command -v julia &>/dev/null
}

# Make juliaup's bin dir stick for future shells (start.sh, new terminals).
persist_juliaup_path() {
    local line='export PATH="$HOME/.juliaup/bin:$PATH"'
    local rc=""
    case "$(basename "${SHELL:-}")" in
        zsh)  rc="${ZDOTDIR:-$HOME}/.zshrc" ;;
        bash) rc="$HOME/.bashrc" ;;
    esac

    echo ""
    echo "NOTE: julia is on PATH via juliaup, but only for this session."
    if [ -n "$rc" ]; then
        if grep -qF "$line" "$rc" 2>/dev/null; then
            echo "      $rc already adds it - run 'source $rc' in open shells."
        else
            printf '\n# Added by MetaManifold install.sh\n%s\n' "$line" >> "$rc"
            echo "      Appended it to $rc - run 'source $rc' or open a new terminal."
        fi
    else
        echo "      Add this to your shell profile so start.sh can find julia:"
        echo "          $line"
    fi
    echo ""
}

if command -v julia &>/dev/null; then
    echo "Found Julia: $(julia --version)"
elif locate_julia; then
    echo "Found Julia (via juliaup, was not on PATH): $(julia --version)"
    persist_juliaup_path
else
    echo "Julia not found. Installing via juliaup..."
    if ! curl -fsSL https://install.julialang.org | sh -s -- --yes; then
        echo ""
        echo "The juliaup installer exited non-zero - usually an existing juliaup"
        echo "install blocking a reinstall. Trying the existing installation..."
    fi

    if locate_julia; then
        echo "Julia available: $(julia --version)"
        persist_juliaup_path
    else
        echo ""
        echo "ERROR: Julia is still not callable after attempting install."
        echo "  Inspect:  juliaup status   ;   ls -la \"$JULIAUP_BIN\""
        echo "  Or install manually: https://julialang.org/downloads/"
        echo "  Then re-run: bash install.sh"
        exit 1
    fi
fi

# Pin this directory to the Julia version the Manifest was resolved with.
JULIA_PIN="$(awk '/^  julia:/{f=1; next} f && /version:/{gsub(/"/, "", $2); print $2; exit}' config/defaults/tool_versions.yml)"
if [ -n "$JULIA_PIN" ] && command -v juliaup &>/dev/null; then
    if juliaup add "$JULIA_PIN" >/dev/null 2>&1 || juliaup status 2>/dev/null | grep -q "$JULIA_PIN"; then
        juliaup override set "$JULIA_PIN" >/dev/null 2>&1 &&
            echo "Pinned Julia $JULIA_PIN for this directory (juliaup override)."
    else
        echo "WARNING: could not install Julia $JULIA_PIN via juliaup; using $(julia --version)."
    fi
fi
echo ""

# R check and politely ask user to do it for us

if command -v Rscript &>/dev/null; then
    echo "Found R:     $(Rscript --version 2>&1 | head -1)"
else
    echo ""
    echo "R is not installed. Please install R for your system and re-run this script."
    echo ""
    if [ "$OS_TYPE" = "Linux" ]; then
        echo "  Ubuntu/Debian:  https://cran.r-project.org/bin/linux/ubuntu/"
        echo "  Fedora/RHEL:    https://cran.r-project.org/bin/linux/fedora/"
        echo "  Quick install:  sudo apt install r-base   (Debian/Ubuntu)"
    elif [ "$OS_TYPE" = "macOS" ]; then
        echo "  macOS pkg:      https://cran.r-project.org/bin/macosx/"
        echo "  Homebrew:       brew install r"
    fi
    echo ""
    exit 1
fi

# Dependency installs by Julia

echo ""
echo "Running install.jl..."
echo ""
julia --project=. install.jl "$@"
