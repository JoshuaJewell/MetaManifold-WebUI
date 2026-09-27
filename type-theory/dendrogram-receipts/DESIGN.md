<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->
<!-- SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell; 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->

# Dendrogram → bridge index → ProofTransport receipt

Status: **exploratory**. The TypeScript and tests are real and pass. The Julia
producer (`export_bridge_index`) is only a proposal (§5).

## 1. The question

> When someone clicks an internal node in the ITOL-style dendrogram, look up
> the bridge index, get the `epistemic_domain` nodes, and show only the
> relevant Agda `ProofTransport` receipt in the right-hand panel, without
> choking the DOM.

## 2. What already exists (so we don't invent it twice)

| Piece | Where | What it gives us |
|---|---|---|
| Taxonomy tree, cumulative rollup, NHX Newick, iTOL bundle | `Protoctist.jl` `src/tree.jl`, `src/io.jl` (`build_taxonomy_tree`, `rollup_counts`, `annotate_newick`, `export_itol_bundle`) | Node ids and a linear post-order rollup, so domain membership can be computed server-side ahead of time |
| Receipts, `epi_status` ladder, standpoints | `EpistemicTypes.jl` (`make_receipt`, `verify_receipt`, `epi_status`) and this repo's `src/core/epistemic.jl` | The `EpiStatus` values and the standpoint (κ) |
| Proof-transport vocabulary | `epistemic-types/.machine_readable/proof-transport/ProofTransport.a2ml` | `statuses`, `modes`, `gaps`, the `Meaning`/`Payload` split |
| Echo fibres | `echo-types`, `EchoTypes.jl`, `proofs/agda/MetaManifold/Evidence/Echo.agda` | Why a clade's rows can be stored as (observed, witnesses) pairs |
| Exact numbers | `src/analysis/exact_summaries.jl`, `proofs/agda/MetaManifold/ExactCounts.agda` | The claims the example receipts are about (`checkedSumOf-is-exact`, …) |
| Clade tree in the UI | `src/analysis/clade_cumulus.jl`, `frontend/src/components/CladeCumulus.tsx` | The node shape the tree view will actually render |

I searched every repo listed above for **`epistemic_domain`**, **"bridge index"** and **JEG**.
None of them appear. This design therefore **defines** the first two (§3). It
treats JEG only as "the host shell", reached through the three-member
`DendrogramHost` interface. If JEG already has its own extension API, only
that adapter has to change.

## 3. The contract (minimal JSON)

There are two documents, fetched at different times.

**Bridge index** (`metamanifold.bridge-index/v1`): one per rendered tree,
loaded with the tree. It is flat maps only, so a click is a hash lookup:

```jsonc
{
  "schema": "metamanifold.bridge-index/v1",
  "tree":     { "newick_sha256": "…", "node_count": 7 },
  "nodes":    { "n2": { "label": "Apicomplexa", "rank": "Class", "domains": ["exacts"] } },
  "domains":  { "exacts": { "kind": "epistemic_domain", "id": "exacts", "label": "…",
                            "standpoint": "…", "epi_status": "Factive",
                            "receipt_ids": ["r-checked-sum-n2", "r-checked-sum-any"] } },
  "receipts": { "r-checked-sum-n2": { "href": "receipts/r-checked-sum-n2.json",
                                      "claim": "CladeTotalIsExact", "status": "Proof",
                                      "artifacts": ["n2"] } }
}
```

- `nodes[id].domains` has **already been rolled up over the clade**. The
  browser never walks descendants.
- Node ids are the NHX `ND=` tags from `annotate_newick`, emitted as
  `data-node-id`.
- `href` must be relative. The parser rejects absolute URLs, which keeps the
  preview-proxy rule and blocks exfiltration.

**Receipt** (`metamanifold.proof-transport-receipt/v1`): one per receipt,
fetched lazily. It carries `holder`, `artifact`, `claim`, `meaning` (the
proposition itself, not just its label), `status`, `mode`, an optional `gap`,
`under` for `ProofUnder`, the `agda` coordinates (module, theorem, source,
`checked_commit`, `safe`), and an optional `exacts` back-link. Full types
are in `frontend/src/type-theory/dendrogram-receipts/types.ts`. Working
examples are in `examples/`.

The parser enforces the a2ml invariants that can be checked at the JSON level:

- A `Proof`/`ProofUnder` cannot carry a gap.
- `OpaqueReceipt` can never be `Proof` (`opaqueNotCertifying`).
- A receipt's `id` must match the href it was requested under. This is the
  `ReplayedArtifact` negative case at the wire layer.

These checks are **not** verification. The panel says so on every card.

## 4. Choosing "the relevant receipt"

This lives in `resolveSelection`, the only place the policy is defined:

1. Collect receipts from every domain on the clicked node, with duplicates
   removed.
2. A receipt is relevant if its `artifacts` names this node, or if it is
   domain-wide (`artifacts: []`). Receipts about *other* nodes are dropped.
3. Sort: node-specific first, then by status strength
   (`Proof > ProofUnder > Receipt > Claimed > Code > Data`), then by id.
4. Render the first one. For the rest, show only a count.

## 5. Not choking the DOM

`controller.ts` follows five rules:

1. **One delegated listener** (`click` and `keydown`) on the tree root, found
   through `closest('[data-node-kind="internal"][data-node-id]')`. There are
   no per-node listeners, so 10k nodes cost the same as 10.
2. **The click handler resolves synchronously** from the in-memory index. It
   does no layout reads and no await before choosing.
3. **Each new selection aborts the previous fetch** (`AbortController`), and a
   generation counter drops any response that still arrives. A late response
   is still cached, because the user often clicks back.
4. **LRU of full receipts** (default 32). Re-clicking the same node does
   nothing.
5. **At most one panel write per animation frame.** The write is a single
   `replaceChildren()` of a detached fragment built with `textContent` only,
   so receipt strings are never parsed as HTML. Moving the selection
   highlight flips one attribute on two elements.

The controller logic (`createReceiptController`) has no DOM dependency, so
`bun test` covers aborts, stale responses and the cache without a browser.

## 6. Julia producer (proposal: Protoctist.jl is the natural home)

Protoctist.jl already builds the tree, does the post-order rollup and writes
the iTOL bundle. The smallest addition would be a function next to
`export_itol_bundle`:

```julia
export_bridge_index(run::ProtistRun, domains, receipts, outdir) -> String
# nodes[ND]    ← union of domain ids over the clade, in the same post-order
#                pass as rollup_counts (linear, not quadratic)
# domains      ← EpistemicTypes standpoint + epi_status per domain
# receipts     ← stubs; the full receipt JSON is written from the Agda build
#                (checked_commit = the commit `just check` passed on)
```

EchoTypes.jl stays test-only, as Protoctist's own design says. It is where a
property test belongs: "every receipt whose artifacts names node *v* is in
the index under *v* or one of *v*'s ancestors".

## 7. Open questions

- **JEG**: what is its extension API? Replace `DendrogramHost` with an
  adapter to it.
- **`epistemic_domain`**: is it meant to be a domain of *claims* (as modelled
  here: exacts, warrant, …) or a taxonomic domain? The shape works for both.
  Only the labels change.
- Should multiple receipts be viewable (a pager), or is "one plus a count"
  the right default?
- Where does the receipt JSON come from in CI: generated from Agda
  `--safe` output in `proofs.yml`, or hand-curated per release?
