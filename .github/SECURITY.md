<!--
SPDX-License-Identifier: AGPL-3.0-only
SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
-->

# Security policy

## Deployment model, stated first

MetaManifold is a **single-user, no-authentication** application: a Julia
server with a browser UI, designed to run on localhost or on a host trusted
by its one operator. The server carries real defences — a browser-origin
allowlist with a DNS-rebinding check, strict Host checking on loopback binds,
inert serving of uploaded report files, and Oxygen's docs and metrics pages
switched off — but **exposing it beyond its intended deployment is the
operator's risk, not a supported configuration.** Denial-of-service against
an internet-exposed instance is out of scope, as is anything that follows
from putting an unauthenticated API on a public network after ignoring the
deployment guidance in the README's *Deployment* section.

## What is in scope

- The HTTP surface: the origin allowlist and Host checks
  (`src/server/server.jl`), the static-file and report-upload serving
  (`_file_middleware`), SPA routing, and every route under
  `src/server/routes/`.
- How configuration and primer values reach the external tools: pipeline
  options and filenames become command lines; injection through them is a
  finding here.
- Path handling: the file routes, report item paths, and study/run names as
  path components.
- Installer and toolchain integrity: every download is SHA256-pinned
  (`config/defaults/tool_versions.yml`); a tampered pin, an unchecked
  substitute, or a way to make the installer fetch unverified content is a
  finding.
- Provenance integrity: `attestation.yml` and the content hashes it records;
  anything that can make a run attest to bytes it did not use.
- The dependency manifests: `Manifest.toml`, `renv.lock`, `frontend/bun.lock`.

## What is out of scope

- Vulnerabilities **in the external tools themselves** (cutadapt, FastQC,
  MultiQC, DADA2, swarm, vsearch, cd-hit-est, MAFFT, trimAl, IQ-TREE, RAxML,
  gappa, and the R packages): report those to their own projects. A way this
  repository invokes them unsafely is in scope; their internal bugs are not.
- Julia, R, bun and DuckDB runtime CVEs: report upstream; this project pins
  versions and will pick up fixes through its normal dependency updates.
- Internet-exposed DoS, as above.

## Supported versions

The default branch only. There are no backported releases; a fix lands on
`main` behind the usual CI gates. The version that matters for any report is
the git commit plus the pinned toolchain — and every run records exactly that
in its `attestation.yml`.

## How to report

**Preferred channel: private vulnerability reporting.** The repository has
it enabled; the *Report a vulnerability* button on the Security tab reaches
the maintainers privately. Use it.

**By email:** contact **developer@joshuajewell.dev** (Jonathan Jewell, the project
owner) with `[security]` in the subject line.

Include, where they apply:

- what you found, and the impact you think it has;
- the endpoint, route, or pipeline stage it affects;
- the deployment mode — localhost, SSH tunnel, or split frontend/API;
- reproduction steps, or the request that triggers it;
- the versions block from the affected run's `attestation.yml` — it names the
  exact commit and toolchain, which usually settles the first two replies.

Reports are acknowledged when triaged, fixed on `main` with the normal CI
gates, and disclosed on the reporter's preferred schedule. Credit is offered
unless you prefer otherwise.

## Mathematical foundations

This project carries machine-checked proofs of some of its statistics, and it
treats **unsound mathematics as a reportable weakness of its own kind**: a
proof that is vacuous, a theorem that does not say what the code relies on, or
a statistical claim the arithmetic does not support. If you find one, that is
a finding — send it to **j.d.a.jewell@open.ac.uk**, where the
proof-theoretic side of the project lives. Questions and open debate about
the foundations belong in the
[mathematical foundations discussion](https://github.com/hyperpolymath/MetaManifold-WebUI/discussions/149);
private concerns about soundness belong in a report.
