// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell; 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
//
// Minimal wire contract for the ITOL-style dendrogram → ProofTransport receipt
// panel. EXPLORATORY: lives in the type-theory work section, not wired into the
// app routes. Design notes: type-theory/dendrogram-receipts/DESIGN.md.
//
// Three payloads, fetched at three different times:
//
//   1. BridgeIndex      — loaded ONCE with the tree. Small. Everything a click
//                         needs to decide *which* receipt to show, O(1).
//   2. ReceiptStub      — inside the index. Enough to rank and label.
//   3. ProofTransportReceipt — fetched LAZILY, one per click, cached.
//
// Nothing here verifies anything. A `status: "Proof"` is a report that the Agda
// check named in `agda` succeeded at `agda.checked_commit`; the browser displays
// that report, it does not re-establish it (epistemic-types §scope:
// "runtime integration" is explicitly not established).

/** Mirrors `[statuses]` in epistemic-types/.machine_readable/proof-transport/ProofTransport.a2ml. */
export const PROOF_STATUSES = ['Data', 'Code', 'Claimed', 'Receipt', 'Proof', 'ProofUnder'] as const
export type ProofStatus = (typeof PROOF_STATUSES)[number]

/** Mirrors `[modes]`. */
export const TRANSPORT_MODES = ['Public', 'Designated', 'IssuerMediated', 'EnvironmentMediated', 'OpaqueReceipt'] as const
export type TransportMode = (typeof TRANSPORT_MODES)[number]

/** Mirrors `[gaps]`. Only the `realized` subset is produced by the Agda verifier today. */
export const GAPS = [
  'TrivialGap', 'DesignatedGap', 'EnvironmentGap', 'IssuerTrustGap', 'OpaqueGap',
  'MissingChecker', 'MissingEvidence', 'MissingContext', 'InvalidEvidence',
] as const
export type Gap = (typeof GAPS)[number]

/**
 * The epi_status ladder shared by EpistemicTypes.jl / src/core/epistemic.jl /
 * Evidence/Warrant.agda (:SansFibre → :Belief → :Warranted → :Factive).
 */
export const EPI_STATUSES = ['SansFibre', 'Belief', 'Warranted', 'Factive'] as const
export type EpiStatus = (typeof EPI_STATUSES)[number]

// ── 1. Bridge index (one JSON document per rendered tree) ────────────────────

export const BRIDGE_INDEX_SCHEMA = 'metamanifold.bridge-index/v1' as const

export interface BridgeIndex {
  schema: typeof BRIDGE_INDEX_SCHEMA
  /** Binds the index to one tree. A click on a tree whose hash differs is refused. */
  tree: { newick_sha256: string; node_count: number }
  /**
   * Internal tree node id → epistemic_domain ids already ROLLED UP over the
   * clade (Protoctist.rollup_counts does the post-order pass server-side), so
   * a click never walks descendants in the browser. Leaves may be omitted.
   * Ids are the NHX `ND=` tag written by Protoctist.annotate_newick.
   */
  nodes: Record<string, BridgeNode>
  domains: Record<string, EpistemicDomain>
  receipts: Record<string, ReceiptStub>
}

export interface BridgeNode {
  /** Taxonomic label, e.g. "Apicomplexa" — the ProofTransport `Artifact`. */
  label: string
  rank: string
  domains: string[]
}

export interface EpistemicDomain {
  kind: 'epistemic_domain'
  id: string
  /** Human label, e.g. "exact-descriptive-summaries". */
  label: string
  /** Standpoint κ (the holder / Agent the claims are indexed by). */
  standpoint: string
  epi_status: EpiStatus
  receipt_ids: string[]
}

export interface ReceiptStub {
  id: string
  /** Relative URL of the full ProofTransportReceipt JSON. Never absolute. */
  href: string
  claim: string
  status: ProofStatus
  /** Tree node ids this receipt's `artifact` is about. Empty = domain-wide. */
  artifacts: string[]
}

// ── 3. Full receipt (fetched per click) ──────────────────────────────────────

export const RECEIPT_SCHEMA = 'metamanifold.proof-transport-receipt/v1' as const

export interface ProofTransportReceipt {
  schema: typeof RECEIPT_SCHEMA
  id: string
  holder: string
  artifact: string
  claim: string
  /** Human rendering of `Meaning holder artifact claim` — the proposition, not the label. */
  meaning: string
  status: ProofStatus
  mode: TransportMode
  /** Present iff verification did not yield Proof/ProofUnder. */
  gap?: Gap
  /** For ProofUnder: the named assumptions the proof is relative to. */
  under?: string[]
  agda: {
    module: string
    theorem: string
    source: string
    checked_commit: string
    safe: boolean
  }
  /** Optional link back to the exact numbers (ExactSummaries) the claim is about. */
  exacts?: { summary_ref: string; policy_fingerprint: string }
}

// ── Host contract ────────────────────────────────────────────────────────────

/**
 * The only three things the extension needs from its host (the JEG shell, or
 * this repo's React app). Kept deliberately tiny so the extension is testable
 * without a browser and portable across hosts.
 */
export interface DendrogramHost {
  /** The element that contains the rendered tree (an <svg> or <div>). */
  treeRoot: Element
  /** The right-hand panel. The controller owns its children exclusively. */
  panel: Element
  fetchJson: (href: string, signal: AbortSignal) => Promise<unknown>
}
