#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# Sourced by CI steps. fetch_pinned downloads a pinned artifact with retries and,
# on failure, annotates the cause. The caller checks the checksum and unpacks.
#
#     source scripts/ci/fetch_pinned.sh
#     fetch_pinned "$VSEARCH_URL" vsearch.tar.gz
#     echo "$VSEARCH_SHA256  vsearch.tar.gz" | sha256sum -c -

# fetch_pinned <url> <output-path>
# Returns 0 on a non-empty download, otherwise curl's exit code.
fetch_pinned() {
    local url="$1" out="$2" rc=0

    # --retry-all-errors covers TLS and receive failures (curl exits 35, 56, 60),
    # which plain --retry skips.
    curl -fsSL --retry 4 --retry-delay 5 --retry-all-errors \
         --connect-timeout 20 --max-time 600 \
         -o "$out" "$url" || rc=$?

    if [[ "$rc" -eq 0 ]]; then
        if [[ -s "$out" ]]; then
            return 0
        fi
        echo >&2 "::error::fetched $url but it is empty"
        return 1
    fi

    local host="${url#https://}"
    host="${host%%/*}"
    case "$host" in
        *:*) ;;
        *)   host="$host:443" ;;
    esac

    case "$rc" in
        35|51|58|60|77|83|90)
            echo >&2 "::error::TLS verification failed (curl exit $rc) fetching $url"
            echo >&2 "::error::Check the host's certificate dates with:"
            echo >&2 "::error::  openssl s_client -connect $host </dev/null 2>/dev/null | openssl x509 -noout -dates -subject"
            echo >&2 "::error::If the artifact has moved, update its url and sha256 in config/defaults/tool_versions.yml."
            ;;
        22)
            echo >&2 "::error::HTTP error status fetching $url (curl exit 22); check the pin in config/defaults/tool_versions.yml"
            ;;
        6|7|28|56)
            echo >&2 "::error::network failure (curl exit $rc) fetching $url after 4 retries"
            ;;
        *)
            echo >&2 "::error::download failed (curl exit $rc) fetching $url"
            ;;
    esac
    return "$rc"
}
