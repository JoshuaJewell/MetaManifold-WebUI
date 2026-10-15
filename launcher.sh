#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# @launcher-deed begin
# ;; SPDX-License-Identifier: MPL-2.0
# (praxis-deed
#   :schema-version  "1.0.0"
#   :canonical-name  "metamanifold-webui-launcher"
#   (artefact :type "launcher" :version "0.1.0"
#             :generator "hand-adapted from launcher-standard 0.6.0's archetype")
#   (app :name "MetaManifold" :display "MetaManifold"
#        :url "https://github.com/hyperpolymath/MetaManifold-WebUI"
#        :archetype "application")
#   (compliance :standard-version "0.6.0"
#               :standards ("launcher-standard.adoc"))
#   (modes :accepted ("--setup" "--start" "--doctor" "--status" "--agent" "--help" "--version"))
#   (platforms :supported ("linux" "macos"))
#   (lifecycle-phases :covered ("provision" "diagnose" "run")
#                     :deferred ("heal")))
# @launcher-deed end
#
# The MetaManifold launcher, in the launcher-standard's mode vocabulary but
# self-contained: this is an application repository, not an RSR-minted one,
# so the modes delegate directly to the repository's own entry points
# (install.sh, start.sh) rather than through provision-set machinery.
#
# The one thing that is this launcher's own is --agent: the git-signing
# agent's persistence is the operator's choice, and every mode of it is
# reversible. The signing key is ~/.ssh/git-signing2; the socket is
# /tmp/git-signing-agent.sock.
#
#   --agent systemd   a systemd user service owns the agent, and the shell
#                      profile exports the socket (a service cannot export
#                      into login shells). Survives reboots; one instance.
#   --agent profile   the shell profile starts the agent if the socket is
#                      dead and exports it. No systemd; works over plain ssh.
#   --agent both      the service plus the spawn-if-dead fallback line.
#   --agent none      removes whatever of the above was installed.
#   --agent status    reports socket, key, and which mechanisms are present.
#
# Each mechanism is written between marked blocks in the profile and a unit
# file, idempotently, with a dated backup of anything edited.

set -uo pipefail

REPO_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
SOCKET="/tmp/git-signing-agent.sock"
KEY="${HOME}/.ssh/git-signing2"
UNIT_NAME="git-signing-agent"
UNIT_FILE="${HOME}/.config/systemd/user/${UNIT_NAME}.service"
PROFILE="${HOME}/.profile"
BLOCK_BEGIN="# >>> metamanifold git-signing agent >>>"
BLOCK_END="# <<< metamanifold git-signing agent <<<"

# Print all arguments as one launcher-prefixed message to standard output.
say()  { printf 'launcher: %s\n' "$*"; }

# Print all arguments as one launcher-prefixed message to standard error.
warn() { printf 'launcher: %s\n' "$*" >&2; }

# Report all arguments through warn and terminate the launcher with status 1.
die()  { warn "$*"; exit 1; }

# Return 0 only when SOCKET is a Unix socket and ssh-add can list identities.
# An agent with no loaded identities also returns non-zero.
agent_alive() {
    [ -S "$SOCKET" ] && SSH_AUTH_SOCK="$SOCKET" ssh-add -l >/dev/null 2>&1
}

# Remove marked agent blocks from PROFILE, following symlinks and preserving
# file permissions. A missing profile succeeds. Return 1 for nested or
# unterminated blocks, or file-operation failures; no backup is made here.
remove_profile_block() {
    [ -f "$PROFILE" ] || return 0
    local target tmp mode
    target="$(readlink -f "$PROFILE")" || return 1
    tmp="$(mktemp "$(dirname "$target")/.profile.tmp.XXXXXX")" || return 1
    if ! awk -v b="$BLOCK_BEGIN" -v e="$BLOCK_END" '
        $0 == b {
            if (skip) exit 1
            skip = 1
            next
        }
        $0 == e {
            if (skip) { skip = 0; next }
            print
            next
        }
        !skip { print }
        END { if (skip) exit 1 }
    ' "$target" > "$tmp"; then
        rm -f "$tmp"
        return 1
    fi
    mode="$(stat -c '%a' "$target" 2>/dev/null || stat -f '%Lp' "$target")" || {
        rm -f "$tmp"
        return 1
    }
    chmod "$mode" "$tmp" || {
        rm -f "$tmp"
        return 1
    }
    mv "$tmp" "$target" || {
        rm -f "$tmp"
        return 1
    }
}

# Back up PROFILE and replace its marked agent block, creating the file if absent.
# "export" sets the socket for future shells; "spawn" also tries to start an
# agent when ssh-add cannot list identities. Neither mode loads the signing key.
# Exit 1 if profile creation, backup or block removal fails; append failures
# return a non-zero status and do not restore the backup.
write_profile_block() {  # $1 = "export" | "spawn"
    [ -f "$PROFILE" ] || touch "$PROFILE" || die "could not create $PROFILE"
    local backup
    backup="$(mktemp "$PROFILE.launcher-backup.$(date +%Y%m%d%H%M%S).XXXXXX")" ||
        die "could not create backup of $PROFILE"
    cp -p "$PROFILE" "$backup" || {
        rm -f "$backup"
        die "backup of $PROFILE failed"
    }
    remove_profile_block || die "could not update $PROFILE"
    {
        printf '%s\n' "$BLOCK_BEGIN"
        if [ "$1" = "spawn" ]; then
            printf 'if [ ! -S "%s" ] || ! SSH_AUTH_SOCK="%s" ssh-add -l >/dev/null 2>&1; then\n' "$SOCKET" "$SOCKET"
            printf '    ssh-agent -a "%s" >/dev/null 2>&1\n' "$SOCKET"
            printf 'fi\n'
        fi
        printf 'export SSH_AUTH_SOCK="%s"\n' "$SOCKET"
        printf '%s\n' "$BLOCK_END"
    } >> "$PROFILE"
}

# Create or overwrite UNIT_FILE for an ssh-agent on SOCKET, without a backup.
# Return the write status; this does not enable or start the service.
write_unit() {
    mkdir -p "$(dirname "$UNIT_FILE")"
    cat > "$UNIT_FILE" <<UNIT
[Unit]
Description=git signing ssh-agent (MetaManifold launcher)

[Service]
ExecStart=/usr/bin/ssh-agent -D -a ${SOCKET}
Restart=on-failure

[Install]
WantedBy=default.target
UNIT
}

# Apply $1: systemd enables the user service and installs a socket export;
# profile installs a fallback that starts an agent; both installs both mechanisms.
# none attempts to disable the service, deletes its unit and removes the profile
# block after backing up the profile. status only reports the observed state.
# Invalid modes, unavailable systemctl when required, service enablement failures
# and profile creation, backup or block removal failures exit 1. In none mode,
# profile backup/removal failures return 1; service disable/reload failures warn.
agent_mode() {
    mode="${1:-}"
    case "$mode" in
        systemd)
            command -v systemctl >/dev/null || die "systemctl not found; use --agent profile"
            write_unit
            systemctl --user daemon-reload
            systemctl --user enable --now "$UNIT_NAME" || die "could not enable ${UNIT_NAME} (try: loginctl enable-linger)"
            write_profile_block export
            say "systemd user service enabled; profile exports SSH_AUTH_SOCK=$SOCKET"
            say "note: over plain ssh, or without a user session, also run: loginctl enable-linger \$USER"
            ;;
        profile)
            write_profile_block spawn
            say "profile now starts the agent if the socket is dead and exports SSH_AUTH_SOCK=$SOCKET"
            ;;
        both)
            command -v systemctl >/dev/null || die "systemctl not found; use --agent profile"
            write_unit
            systemctl --user daemon-reload
            systemctl --user enable --now "$UNIT_NAME" || die "could not enable ${UNIT_NAME}"
            write_profile_block spawn
            say "systemd service plus spawn-if-dead profile fallback installed"
            ;;
        none)
            if command -v systemctl >/dev/null 2>&1; then
                systemctl --user disable --now "$UNIT_NAME" >/dev/null 2>&1 ||
                    warn "could not disable $UNIT_NAME"
            else
                warn "systemctl not found; user unit was not disabled"
            fi
            rm -f "$UNIT_FILE"
            if command -v systemctl >/dev/null 2>&1; then
                systemctl --user daemon-reload ||
                    warn "could not reload the user manager"
            fi
            if [ -f "$PROFILE" ]; then
                backup="$(mktemp "$PROFILE.launcher-backup.$(date +%Y%m%d%H%M%S).XXXXXX")" || {
                    warn "could not create backup of $PROFILE; profile block was not removed"
                    return 1
                }
                if ! cp -p "$PROFILE" "$backup"; then
                    rm -f "$backup"
                    warn "could not back up $PROFILE; profile block was not removed"
                    return 1
                fi
            fi
            if ! remove_profile_block; then
                warn "could not remove the profile block"
                return 1
            fi
            say "removed unit and profile block; any agent already running is left alone"
            ;;
        status)
            agent_alive && say "socket alive: $SOCKET" || say "socket dead or absent: $SOCKET"
            [ -f "$KEY" ] && say "key present: $KEY" || say "key MISSING: $KEY"
            [ -f "$UNIT_FILE" ] && say "systemd unit: installed" || say "systemd unit: not installed"
            grep -q "$BLOCK_BEGIN" "$PROFILE" 2>/dev/null && say "profile block: present" || say "profile block: absent"
            ;;
        *)
            die "usage: launcher.sh --agent systemd|profile|both|none|status"
            ;;
    esac
}

# Report toolchain and signing readiness, then exit 1 if julia, bun, R or git
# is missing from PATH, or agent_alive fails; otherwise exit 0. Missing tool
# configuration or a signing-key check failure only produces a warning.
doctor() {
    fail=0
    for tool in julia bun R git; do
        if command -v "$tool" >/dev/null; then
            say "ok: $tool on PATH"
        else
            warn "MISSING: $tool not on PATH"
            fail=1
        fi
    done
    [ -f "$REPO_DIR/config/tools.yml" ] && say "ok: config/tools.yml exists" \
        || warn "config/tools.yml missing (run: ./launcher.sh --setup)"
    if agent_alive; then
        say "ok: signing agent alive at $SOCKET"
        SSH_AUTH_SOCK="$SOCKET" ssh-add -l 2>/dev/null | grep -q "$(ssh-keygen -lf "${KEY}.pub" 2>/dev/null | awk '{print $2}')" \
            && say "ok: signing key loaded" || warn "signing key not loaded (run: ssh-add $KEY)"
    else
        warn "signing agent not reachable at $SOCKET (see: ./launcher.sh --agent status)"
        fail=1
    fi
    exit "$fail"
}

# Print the launcher header and available modes to standard output.
help() {
    sed -n '2,40p' "${BASH_SOURCE[0]}" | grep '^#' | sed 's/^# \{0,1\}//'
    cat <<'MODES'
Modes:
  --setup            install everything (delegates to install.sh; args pass through)
  --start            serve the web UI (delegates to start.sh)
  --doctor           check the toolchain and the signing agent
  --status           whether the server looks reachable on its port
  --agent MODE       systemd | profile | both | none | status  (the signing agent)
  --help             this text
  --version          the launcher's version
MODES
}

# Print the launcher and launcher-standard versions to standard output.
version() { echo "metamanifold-webui-launcher 0.1.0 (launcher-standard 0.6.0, hand-adapted)"; }

# Dispatch the first argument as a launcher mode. Setup replaces this process
# with install.sh and forwards remaining arguments; start executes start.sh.
# Status checks the configured local port (default 8080) with a two-second
# timeout; HTTP error responses count as reachable and connection failure does
# not set a failure exit status. Missing or unknown modes print help and exit 1.
main() {
    case "${1:-}" in
        --setup)  shift; exec bash "$REPO_DIR/install.sh" "$@" ;;
        --start)  exec bash "$REPO_DIR/start.sh" ;;
        --doctor) doctor ;;
        --status)
            status_url="http://127.0.0.1:${JULIA_METAMANIFOLD_PORT:-8080}/"
            if curl -s -o /dev/null --max-time 2 "$status_url"; then
                say "something is answering on $status_url"
            else
                say "nothing answering on $status_url (start with: ./launcher.sh --start)"
            fi
            ;;
        --agent)  shift; agent_mode "${1:-}" ;;
        --help|-h) help ;;
        --version) version ;;
        *)        help; exit 1 ;;
    esac
}

main "$@"
