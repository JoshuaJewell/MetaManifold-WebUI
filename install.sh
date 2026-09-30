#!/usr/bin/env bash
#
# Bootstrap installer for MetabarcodingPipeline
#
# Checks for Julia and R, installs the pinned Julia release if Julia is missing,
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

# The pinned Julia distribution, unpacked by install_pinned_julia. start.sh
# looks here too.
PINNED_JULIA_DIR="bin/julia"
PINS_FILE="config/defaults/tool_versions.yml"

# Print one field ("url" or "sha256") of runtimes.julia.archives.<platform>.
julia_archive_field() {
    awk -v plat="$1" -v key="$2" '
        /^runtimes:/            { r = 1; next }
        r && /^[^ #]/           { r = 0 }
        r && /^  julia:/        { j = 1; next }
        j && /^  [^ ]/          { j = 0 }
        j && $1 == plat ":"     { p = 1; next }
        p && /^      [^ ]/      { p = 0 }
        p && $1 == key ":"      { gsub(/"/, "", $2); print $2; exit }
    ' "$PINS_FILE"
}

# Download the Julia release pinned in tool_versions.yml, check it against the
# recorded sha256 before unpacking anything, and unpack it into
# $PINNED_JULIA_DIR. This replaces piping a live installer into sh, which ran
# whatever that URL served on the day with no checksum at all.
install_pinned_julia() {
    local arch platform url sha tmp actual
    arch="$(uname -m)"
    case "$arch" in
        x86_64|amd64)  arch="x86_64"  ;;
        aarch64|arm64) arch="aarch64" ;;
    esac
    case "$OS_TYPE" in
        Linux) platform="linux-$arch" ;;
        macOS) platform="macos-$arch" ;;
    esac

    url="$(julia_archive_field "$platform" url)"
    sha="$(julia_archive_field "$platform" sha256)"
    if [ -z "$url" ] || [ -z "$sha" ]; then
        echo "ERROR: no pinned Julia archive for $platform in $PINS_FILE."
        return 1
    fi

    mkdir -p bin
    tmp="$(mktemp -d bin/.julia-download.XXXXXX)"
    echo "Downloading $url"
    if ! curl -fsSL --proto '=https' --proto-redir '=https' --retry 3 \
              -o "$tmp/julia.tar.gz" "$url"; then
        rm -rf "$tmp"
        echo "ERROR: download failed: $url"
        return 1
    fi

    if command -v sha256sum &>/dev/null; then
        actual="$(sha256sum "$tmp/julia.tar.gz" | awk '{print $1}')"
    else
        actual="$(shasum -a 256 "$tmp/julia.tar.gz" | awk '{print $1}')"
    fi
    if [ "$actual" != "$sha" ]; then
        rm -rf "$tmp"
        echo "ERROR: checksum mismatch for $url"
        echo "  expected $sha"
        echo "  got      $actual"
        return 1
    fi

    mkdir "$tmp/julia"
    tar -xzf "$tmp/julia.tar.gz" -C "$tmp/julia" --strip-components=1
    if [ ! -x "$tmp/julia/bin/julia" ]; then
        rm -rf "$tmp"
        echo "ERROR: the archive did not contain bin/julia."
        return 1
    fi
    rm -rf "$PINNED_JULIA_DIR"
    mv "$tmp/julia" "$PINNED_JULIA_DIR"
    rm -rf "$tmp"
}

# Precedence: a pinned release already in bin/julia, then julia on PATH, then
# juliaup, and only then a fresh pinned install. Whichever is chosen must be
# the pinned version; the check after the juliaup override below enforces that.
if [ -x "$PINNED_JULIA_DIR/bin/julia" ]; then
    PATH="$PWD/$PINNED_JULIA_DIR/bin:$PATH"
    echo "Found Julia (pinned, in $PINNED_JULIA_DIR): $(julia --version)"
elif command -v julia &>/dev/null; then
    echo "Found Julia: $(julia --version)"
elif locate_julia; then
    echo "Found Julia (via juliaup, was not on PATH): $(julia --version)"
    persist_juliaup_path
else
    echo "Julia not found. Installing the pinned release into $PINNED_JULIA_DIR..."
    if install_pinned_julia; then
        PATH="$PWD/$PINNED_JULIA_DIR/bin:$PATH"
        echo "Julia available: $(julia --version)"
        echo "start.sh finds it in $PINNED_JULIA_DIR; no shell profile change is needed."
    else
        echo ""
        echo "ERROR: Julia could not be installed."
        echo "  Install it manually (https://julialang.org/downloads/), or install"
        echo "  juliaup, then re-run: bash install.sh"
        exit 1
    fi
fi

# Pin this directory to the Julia version the Manifest was resolved with.
JULIA_PIN="$(awk '/^  julia:/{f=1; next} f && /version:/{gsub(/"/, "", $2); print $2; exit}' "$PINS_FILE")"
if [ -n "$JULIA_PIN" ] && command -v juliaup &>/dev/null; then
    if juliaup add "$JULIA_PIN" >/dev/null 2>&1 || juliaup status 2>/dev/null | grep -q "$JULIA_PIN"; then
        juliaup override set "$JULIA_PIN" >/dev/null 2>&1 &&
            echo "Pinned Julia $JULIA_PIN for this directory (juliaup override)."
    else
        echo "WARNING: could not install Julia $JULIA_PIN via juliaup."
    fi
fi

# Whatever was found, Manifest.toml was resolved with $JULIA_PIN. A julia of any
# other version (a distro package, a stale bin/julia after the pin moved, a
# juliaup that could not be overridden) is replaced here by the pinned release.
JULIA_FOUND="$(julia --version 2>/dev/null | awk '{print $3}')"
if [ -n "$JULIA_PIN" ] && [ "$JULIA_FOUND" != "$JULIA_PIN" ]; then
    echo "Julia ${JULIA_FOUND:-(unknown)} does not match the pinned $JULIA_PIN."
    echo "Installing the pinned release into $PINNED_JULIA_DIR..."
    if install_pinned_julia; then
        PATH="$PWD/$PINNED_JULIA_DIR/bin:$PATH"
        echo "Using $(julia --version) from $PINNED_JULIA_DIR."
    else
        echo "WARNING: could not install Julia $JULIA_PIN; continuing with $JULIA_FOUND."
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
