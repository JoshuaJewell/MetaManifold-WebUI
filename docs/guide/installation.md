<!-- SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# Installation and deployment

## Prerequisites

- **Julia** 1.12 (installed by `install.sh` if missing, and pinned to the Manifest version for this directory)
- **R** >= 4.0 (required for the DADA2 stage and NMDS/PERMANOVA analysis)
  - Ubuntu/Debian: `sudo apt install r-base`
  - macOS: `brew install r` or [CRAN package](https://cran.r-project.org/bin/macosx/)
- **bun** builds the frontend. The installer uses a bun on PATH at the version pinned in `config/defaults/tool_versions.yml`, or downloads that release into `bin/`.

## Installation

```bash
git clone https://github.com/JoshuaJewell/MetaManifold-WebUI.git
cd MetaManifold-WebUI
bash install.sh
```

`install.sh` installs Julia (via juliaup) if missing, checks for R, then hands
off to `install.jl`, which installs the Julia dependencies, locates or downloads
each external tool (cutadapt, FastQC, MultiQC, vsearch, cd-hit-est, and the
phylogeny tools MAFFT, trimAl, IQ-TREE, RAxML and gappa), and reproduces the R
environment from `renv.lock`. It ends with a summary listing every dependency as
`OK`, `SKIPPED`, or `ACTION NEEDED`. Tool paths can be set in
`config/tools.yml`. MAFFT, IQ-TREE and RAxML can be skipped on a machine that
sends those stages to a server (see
[Threads and remote execution](configuration.md#threads-and-remote-execution-top-level-of-pipelineyml)).

Installing R and the system `-dev` headers the R packages compile against needs
root. If R is missing, the installer stops and prints how to install it. If
headers are missing, the summary prints the command that installs them; run it
and rerun `install.sh`.

To update:
```bash
bash install.sh --update
```

### Reproducing the R environment

The R-side dependencies (DADA2, vegan, and their transitive packages) are
pinned with [`renv`](https://rstudio.github.io/renv/); `renv.lock` is the source
of truth and the committed `.Rprofile` activates the project-local library on
any `R`/`Rscript` invocation from the repository root. `install.sh` runs the
restore; to run it by hand:

```bash
Rscript -e 'renv::restore(prompt = FALSE)'
```

To add or upgrade a package, install it inside the project
(`renv::install(...)`) and record the change with `renv::snapshot()`.

### Tool paths

`install.sh` generates `config/tools.yml`, one local path per tool:

```yaml
vsearch:
  path: "/home/user/software/vsearch"
```

`config/defaults/tools.yml` contains the full config format. Stages that run on a
server are configured in the `remote` block of `pipeline.yml`.

## Running the server

| Environment variable      | Default           | Description   |
| ------------------------- | ----------------- | ------------- |
| `JULIA_METAMANIFOLD_PORT` | `8080`            | Server port   |
| `JULIA_METAMANIFOLD_ROOT` | working directory | Project root  |
| `JULIA_THREADS`           | `8`               | Julia threads |
| `JULIA_METAMANIFOLD_ALLOWED_ORIGINS` | none | Extra browser origins allowed to call the API, comma-separated (the server's own origin and the Vite dev server at `localhost:5173` are always allowed) |

Or run the Julia server directly:

```bash
julia --project=. scripts/serve.jl
```

The server is part of the `MetaManifold` package, so the first start after an
install or a source change compiles the package. That compile runs a workload of
the requests a first page visit makes, which takes a few minutes and caches
their native code. While developing, skip the workload with a
`LocalPreferences.toml` (git-ignored) in the repository root:

```toml
[MetaManifold]
precompile_workload = false
```

The server then makes those requests itself in the background as it starts, so the first page visit does not wait on compilation.

## Deployment

### Local (single machine)

```bash
bash start.sh
```

Open `http://localhost:8080`. The backend serves the frontend.

### Split deployment (frontend and API on different origins)

To serve the frontend from a different address than the API, set `apiBase` in
`frontend/public/config.json` to the API's address, and list the frontend's
address when starting the server:

```bash
JULIA_METAMANIFOLD_ALLOWED_ORIGINS=https://frontend.example bash start.sh
```

The server accepts browser requests from its own origin (`http://localhost:<port>`
and `http://127.0.0.1:<port>`), from the Vite dev server (`localhost:5173`), and
from the origins listed in `JULIA_METAMANIFOLD_ALLOWED_ORIGINS`, separated by
commas. An origin has to match exactly, scheme and port included:
`https://frontend.example` does not cover `http://frontend.example` or
`https://frontend.example:8443`. A page on any other origin gets
`403 Forbidden: origin not allowed` on every request that carries an `Origin`
header, which is every write and every cross-origin read.

This also applies to an SSH tunnel on a different local port (e.g.
`http://localhost:9000`): the browser's origin is then `http://localhost:9000`,
not the server's `http://localhost:8080`, so start the server with
`JULIA_METAMANIFOLD_ALLOWED_ORIGINS=http://localhost:9000`.
