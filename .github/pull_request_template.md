## Summary

<!-- What does this PR do, and why? One logical change per PR. -->

Closes #

## Type of change

- [ ] 🐛 Bug fix (non-breaking change that fixes an issue)
- [ ] ✨ New feature (non-breaking change that adds functionality)
- [ ] 💥 Breaking change (would change existing behaviour, a config key or an API response)
- [ ] 📊 Statistics fix (a method returned a wrong number, or a result it could not support)
- [ ] 📖 Documentation
- [ ] 🧹 Refactor / tech debt (behaviour-preserving)
- [ ] ⚡ Performance
- [ ] 🔧 Build / CI / tooling

## How has this been verified?

<!-- Name the commands you ran and what they reported. For example:
     julia --project=. test/runtests.jl
     julia --project=. -t 2 test/runtests.jl --integration --server
     cd frontend && bun run test -->

Tested with Julia `…`, R `…` and bun `…`.

## What this leaves out

<!-- Anything related that this PR does not do, and why. -->

## Checklist

- [ ] My branch starts at the current tip of `main` (`behind_by` is `0`; see CONTRIBUTING.md) and has no merge commits from `main`.
- [ ] I ran `julia --project=. test/runtests.jl` (and `cd frontend && bun run test` for frontend changes) and they pass.
- [ ] CI ran on this PR and passes.
- [ ] No generated or vendored files, and no lockfile change the work does not need (or the reason is given above).
- [ ] New code is AGPL-3.0, or MPL-2.0 with an SPDX identifier in each file I wrote. I left existing file headers, `LICENSE`, `CITATION.cff` and the licence and acknowledgement sections of `README.md` unchanged.
- [ ] Every number shown to a user is computed from their data. When a method cannot run it raises an explicit error.
- [ ] Docs are updated, and no claim in them overstates what the code does.

## Notes for reviewers

<!-- Anything that needs special attention, follow-up, or context. To be considered for authorship in CITATION.cff, say so here. -->
