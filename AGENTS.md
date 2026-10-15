# AGENTS.md

Instructions for coding agents working in this repository. Read the contribution
rules in `.github/CONTRIBUTING.md` first; this file records only the traps that
are invisible until they bite.

## renv.lock and CI's hardcoded R-package lists are coupled

CI installs R from `renv.lock` (`renv::restore()`), but `.github/workflows/ci.yml`
names the expected packages in TWO hardcoded verification lists:

- the `Install R packages` step's `require()` check ("failed to install"), and
- the `Verify R packages visible from RCall` step.

Removing a package from `renv.lock` therefore fails the run at install
verification, nowhere near the code that stopped using it. MEASURED on
JoshuaJewell/MetaManifold-WebUI#57 (2026-10-15): dplyr left the lock and both
checks failed until the lists were updated in the same PR. Whenever a change
touches `renv.lock`, grep `.github/workflows/ci.yml` for every added or removed
package name and update both lists in the same commit.

Corollary, measured the same way: `tibble` cannot leave the lock while dada2 is
in it. `dada2 -> ggplot2 -> tibble` is dada2's own dependency chain, independent
of anything this project's code does; the lock's `Requirements` arrays are not
evidence (only a handful of records carry one).

## Commit signing needs the SSH agent, not a passphrase prompt

`commit.gpgsign` is on, with a passphrase-protected SSH key. A non-interactive
shell has no askpass and no terminal for a prompt, so signing fails with a
misleading "failed to write commit object" unless the agent is reachable:
`SSH_AUTH_SOCK=/tmp/lcb-ssh-agent.sock` (the fixed socket; `~/.profile` exports
it when the socket exists). Never work around this with
`commit.gpgsign=false` — check the agent first with
`SSH_AUTH_SOCK=/tmp/lcb-ssh-agent.sock ssh-add -l`.
