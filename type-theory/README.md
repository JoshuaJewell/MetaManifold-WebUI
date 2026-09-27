<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell; 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->

# Type-theory work section (exploratory)

A space for trying out ideas where the estate's type-theory work (echo-types,
epistemic-types, `proofs/agda/`) meets the MetaManifold UI. **Nothing here
is wired into the app routes, the server, or the release.** When something here
is ready, it moves into `src/`, `frontend/src/views/` or `proofs/` through a
normal PR. Until then, treat it as an example of what this work could be used
for, not a promise.

Layout:

| Path | What |
|---|---|
| `type-theory/<topic>/DESIGN.md` | The argument, the contract, open questions |
| `type-theory/<topic>/examples/` | JSON fixtures, used by the tests |
| `frontend/src/type-theory/<topic>/` | TypeScript (checked by the normal `tsc` gate) |
| `frontend/tests/unit/type-theory-*.test.ts` | `bun test` coverage |

Topics:

- [`dendrogram-receipts/`](dendrogram-receipts/DESIGN.md): clicking an internal
  node of an ITOL-style dendrogram looks up the bridge index, gets the
  `epistemic_domain`s for that clade, and shows one Agda `ProofTransport`
  receipt in the right-hand panel.
