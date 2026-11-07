<!-- SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell -->
<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# Development

## REST API

The server exposes a REST API under `/api/v1/`. Key endpoint groups:

| Group          | Endpoints                                                                                                     | Description                                                                              |
| ----------------| ---------------------------------------------------------------------------------------------------------------| ------------------------------------------------------------------------------------------|
| Studies        | `GET/POST/DELETE /studies`, `POST .../rename`                                                                 | List, create, rename, delete studies                                                     |
| Groups         | `POST/DELETE /studies/{study}/groups`, `POST .../rename`                                                      | Create, rename, delete groups                                                            |
| Runs           | `GET/POST/DELETE /studies/{study}/runs`, `POST .../rename`                                                    | List, create, rename, delete runs                                                        |
| Config         | `GET/PATCH/DELETE .../config`, `GET .../config/overrides`                                                     | Read and edit config at any cascade level; list downstream overrides                     |
| Primers        | `GET /primers`, `GET /primers/document`, `PUT /primers`                                                       | List pair names; read and replace the whole primers document (validated before writing)  |
| Pipeline       | `POST .../pipeline`, `POST .../stages/{stage}`                                                                | Launch full-study, single-run, or individual-stage jobs                                  |
| Jobs           | `GET /jobs`, `GET/DELETE /jobs/{id}`                                                                          | Monitor and cancel running pipeline jobs                                                 |
| Events         | `GET /events`                                                                                                 | Server-sent event stream of real-time job and stage updates                              |
| Results        | `GET/POST/DELETE .../results/tables/...`                                                                      | List, query, filter, save, export (`.xlsx`), and delete tables; OTU member drill-down    |
| QC             | `GET .../results/qc`, `GET .../results/dada2`                                                                 | MultiQC report metadata and DADA2 figures, logs, stats                                   |
| Analysis       | `POST .../runs/{run}/analysis/{alpha,chart}`, `GET .../runs/{run}/analysis/{pipeline-stats,ranks}`            | Per-run charts and rank discovery                                                        |
| Cross-run      | `POST /studies/{study}/analysis/{alpha,chart,chart-facet,nmds,permanova,venn,publication-tables}`             | Comparison, faceted charts, NMDS, PERMANOVA, taxon overlap and publication tables across runs |
| Composition    | `POST .../runs/{run}/composition/summary`, `POST .../composition/{source}/query`, `POST .../composition/{source}/distinct/{column}` | Organism-composition summaries and tables for a run                             |
| Composition library | `GET /composition`, `POST/DELETE /composition/{filters,sets}/{name}`, `GET /category-sets`, `POST/DELETE /category-sets/{name}` | Edit the filters and category sets in `config/composition.yml`               |
| Reference trees | `GET/POST /reference-trees`, `GET/PUT/DELETE /reference-trees/{id}`, `GET/PUT .../fasta`, `POST .../run`, `GET .../log/{step}`, `GET .../qc/{align,trim}`, `GET .../alignment/{raw,trimmed}`, `GET .../treefile`, `POST .../trim-preview` | Build and inspect the shared reference trees |
| Placements     | `GET/POST /studies/{study}/placements`, `GET/PUT/DELETE .../placements/{id}`, `POST .../placements/preview`, `GET/PUT .../fasta/queries`, `POST .../run`, `GET .../log/{step}`, `GET .../qc/{align,trim}`, `GET .../alignment/{raw,trimmed}`, `POST .../trim-preview` | Place a study's sequences on a reference tree |
| Trees          | `GET/POST /studies/{study}/trees`, `GET/DELETE /studies/{study}/trees/{file}`, `PUT /studies/{study}/trees/{file}/view` | Upload, read and delete Newick and jplace files; save each tree's view state          |
| Filter presets | `GET/POST/DELETE /filter-presets`                                                                             | Save, list, delete reusable table filters                                                |
| Databases      | `GET /databases`, `GET /databases/document`, `PUT /databases`, `POST /databases/{key}/download`               | List and download taxonomy databases; read and replace the whole databases document (validated before writing, returns advisory warnings) |
| System         | `POST /init`, `GET /capabilities`                                                                             | Initialise project directories; report server capabilities (e.g. R availability)         |

All responses are JSON. Analysis endpoints return Plotly chart specifications.

## Running tests

```bash
julia --project=. test/runtests.jl                  # unit tests
julia --project=. test/runtests.jl --integration    # adds the MiSeq SOP pipeline run
julia --project=. test/runtests.jl --server         # adds the HTTP server smoke tests
```

The integration run assigns taxonomy against PR2 by default, which needs about 16 GiB of memory, because DADA2 holds about 262 KB per distinct lineage. Set `CI_TAXONOMY_DB=<fasta>` to use a smaller reference (CI cuts one with `config/ci/subset_pr2.jl`), or `CI_SKIP_TAXONOMY=1` to skip the stage. The full suite loads R and every package and needs several GB of memory; on a small machine run a few files at a time. `cd frontend && bun run test` runs the frontend unit tests.

Name unit test files to run only those: `julia --project=. test/runtests.jl routes trees`. `dev/dev.sh` runs the development server, rebuilds and tests one job at a time, stopping and restarting the server around builds and tests (`dev/dev.sh` alone lists its commands). Deferred work is listed in `dev/TODO.md`.

## Architecture

```
frontend/           TypeScript + React + Vite (SPA)
src/
  core/             Types, config cascade, validation, DuckDB store, logging
  pipeline/         Pipeline stages (cutadapt, dada2, swarm, vsearch, cd-hit-est, merge_taxa, phylogeny)
  analysis/         Diversity metrics + Plotly chart builders
  server/           Oxygen.jl HTTP server (MetaManifold.Server) and its precompile workload
    routes/         REST API route handlers
scripts/            serve.jl (server entry point) and one-off maintenance scripts
config/             Default configs, filters, CI fixtures
data/               Input FASTQs (user-managed)
projects/           Pipeline outputs (generated)
```

Each pipeline stage returns a typed result (`TrimmedReads`, `ASVResult`, `OTUResult`, `TaxonomyHits`, `MergedTables`) and is skipped when its outputs are up to date (mtime-based for files, content-hash-based for configuration). Rerunning after a config change only re-executes the minimum necessary stages.
