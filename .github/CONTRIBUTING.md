# Contributing to MetaManifold

This file sets out what a pull request must look like, how licensing and credit
work, and which design decisions are settled. Pull requests that miss the
requirements under [Pull requests](#pull-requests) are closed without review.

## Issues

Issues track defects and features in this repository. Raise issues about a fork's
own code on that fork. Search the existing issues before opening one.

A bug report needs the command or steps you ran, what you expected, what
happened, and the Julia, R and browser versions involved.

To change an architectural decision, open an issue describing the change before
writing any code.

## Pull requests

Each requirement is checked before the code is read.

1. **Branch from the current `main`.** The merge base of your branch and `main`
   must be the tip of `main`. Check before pushing:

   ```
   gh api repos/JoshuaJewell/MetaManifold-WebUI/compare/main...<head-sha> \
     --jq '{behind_by, merge_base: .merge_base_commit.sha}'
   ```

   `behind_by` must be `0`. Rebase when it is not.
2. **No merge commits from `main`** and no automated conflict-resolution commits.
   A branch that cannot be re-applied cleanly on top of `main` is not ready.
3. **One logical change per PR.**
4. **CI must run and pass.** A PR with no reported checks is not reviewed.
5. **Write the description:** what the change does, the issues it closes, how to
   run its tests, what it leaves out, and the Julia, R and bun versions you tested
   with.
6. **No generated or vendored files:** build output, `node_modules`, downloaded
   tool archives, or lockfile changes the work does not need. When a lockfile
   change is needed, the description says why.
7. **Once review has started, add commits,** and rebase only when asked to.
   Squashing happens on merge.

## Licence and copyright

Source code is licensed under the AGPL-3.0 (`LICENSE`). `README.md` is licensed
under CC BY-SA 4.0, as its Licence section states. A contribution is offered
under the same terms, with one option: code you write yourself may be licensed
under MPL-2.0. Mark each such file with an SPDX identifier and leave out the
"Incompatible With Secondary Licenses" notice, since that notice prevents
combining the file with AGPL code. Checking that your licence choice combines
correctly is your responsibility. Documentation you write is CC BY-SA 4.0.

In a contribution PR:

- leave `LICENSE`, `CITATION.cff` and the licence and acknowledgement sections of
  `README.md` unchanged, and open an issue if one of them looks wrong;
- add no `NOTICE` file and no `LICENSES/` directory, and do not restate this
  project's licensing in your own files;
- keep existing file headers exactly as they are, and add your copyright line only
  to files you wrote;
- put third-party attributions, such as the DADA2 tutorial used under CC BY 4.0,
  in the acknowledgements section of `README.md`.

## Citation and authorship

The maintainer keeps `CITATION.cff`. Being added to its `authors` is an
authorship decision made by the maintainer, and owning copyright in a
contributed file does not confer it. To be considered for
authorship, say so in the PR description. Your copyright in your contributions is
unaffected by whether you are listed. Credit short of authorship goes in the
acknowledgements section of `README.md`.

## Settled decisions

A pull request does not reopen these. An issue can discuss them.

- **HTTP layer: Oxygen.jl,** pinned in `Manifest.toml`, with the server in
  `src/server/`. Genie, Stipple and other Julia web frameworks are out of scope,
  as is any second Julia web environment.
- **Frontend: React 18 and TypeScript,** built by Vite and managed with bun, in
  `frontend/`. There is one UI stack.
- **Statistics: R through RCall.** Every R package used is pinned in `renv.lock`.
  When a method needs a package missing from the lockfile, the code raises a clear
  error.
- **Storage: DuckDB,** through `src/core/duckdb_store.jl`.

## Statistics integrity

- Every number shown to a user is computed from their data. Placeholder and mock
  values have no place in an analysis path.
- When a method cannot run (a missing package, invalid input, a model that does
  not converge), raise an explicit error. Substituting a cheaper method, or
  returning an empty result that reads as a negative finding, is a defect.
- Tests compare computed values with known references.

## Toolchain and pins

- `config/defaults/tool_versions.yml` pins Julia, R, bun and the external tools,
  and `install.jl` and CI both install from it. The Julia version is repeated in
  the CI workflow. Change a pin by editing that file, and expect the change to be
  discussed.
- R packages are pinned in `renv.lock`.
- CI job names carry no version numbers, so that a pin change leaves required
  status checks intact.

Run the tests as CI does:

```
julia --project=. -t 2 test/runtests.jl --integration --server
cd frontend && bun run test
```

## Contributions written with coding agents

The same rules apply, and a named person answers for every line of the diff.

## Getting help

Ask in an issue, or comment on an existing one.
