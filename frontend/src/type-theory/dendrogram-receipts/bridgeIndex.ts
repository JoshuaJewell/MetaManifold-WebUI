// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell; 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
//
// Pure functions: parse the two JSON payloads at the boundary, and decide which
// single receipt a click on an internal node should show. No DOM, no fetch —
// so the whole selection policy is unit-testable under `bun test`.

import {
  BRIDGE_INDEX_SCHEMA, EPI_STATUSES, GAPS, PROOF_STATUSES, RECEIPT_SCHEMA, TRANSPORT_MODES,
  type BridgeIndex, type EpistemicDomain, type ProofStatus, type ProofTransportReceipt, type ReceiptStub,
} from './types'

export class BridgeContractError extends Error {
  override name = 'BridgeContractError'
}

// ── boundary parsing (unknown → typed, or throw) ────────────────────────────

type Obj = Record<string, unknown>
const isObj = (v: unknown): v is Obj => typeof v === 'object' && v !== null && !Array.isArray(v)
const str = (o: Obj, k: string, at: string): string => {
  const v = o[k]
  if (typeof v !== 'string') throw new BridgeContractError(`${at}.${k}: expected string`)
  return v
}
const strArr = (o: Obj, k: string, at: string): string[] => {
  const v = o[k]
  if (!Array.isArray(v) || !v.every((x) => typeof x === 'string')) {
    throw new BridgeContractError(`${at}.${k}: expected string[]`)
  }
  return v as string[]
}
const oneOf = <T extends string>(allowed: readonly T[], v: unknown, at: string): T => {
  if (typeof v === 'string' && (allowed as readonly string[]).includes(v)) return v as T
  throw new BridgeContractError(`${at}: ${JSON.stringify(v)} not in {${allowed.join(', ')}}`)
}
const record = (o: Obj, k: string, at: string): Obj => {
  const v = o[k]
  if (!isObj(v)) throw new BridgeContractError(`${at}.${k}: expected object`)
  return v
}

/** Receipt hrefs must stay relative (the preview/proxy rule, and no exfiltration). */
function relativeHref(h: string, at: string): string {
  if (/^[a-z][a-z0-9+.-]*:/i.test(h) || h.startsWith('//')) {
    throw new BridgeContractError(`${at}: href must be relative, got ${h}`)
  }
  return h
}

export function parseBridgeIndex(raw: unknown): BridgeIndex {
  if (!isObj(raw)) throw new BridgeContractError('index: expected object')
  if (raw['schema'] !== BRIDGE_INDEX_SCHEMA) {
    throw new BridgeContractError(`index.schema: expected ${BRIDGE_INDEX_SCHEMA}`)
  }
  const tree = record(raw, 'tree', 'index')
  const nodeCount = tree['node_count']
  if (typeof nodeCount !== 'number' || !Number.isInteger(nodeCount) || nodeCount < 0) {
    throw new BridgeContractError('index.tree.node_count: expected non-negative integer')
  }

  const receipts: Record<string, ReceiptStub> = {}
  for (const [id, r] of Object.entries(record(raw, 'receipts', 'index'))) {
    const at = `receipts.${id}`
    if (!isObj(r)) throw new BridgeContractError(`${at}: expected object`)
    receipts[id] = {
      id,
      href: relativeHref(str(r, 'href', at), `${at}.href`),
      claim: str(r, 'claim', at),
      status: oneOf(PROOF_STATUSES, r['status'], `${at}.status`),
      artifacts: strArr(r, 'artifacts', at),
    }
  }

  const domains: Record<string, EpistemicDomain> = {}
  for (const [id, d] of Object.entries(record(raw, 'domains', 'index'))) {
    const at = `domains.${id}`
    if (!isObj(d)) throw new BridgeContractError(`${at}: expected object`)
    if (d['kind'] !== 'epistemic_domain') throw new BridgeContractError(`${at}.kind: expected epistemic_domain`)
    const receiptIds = strArr(d, 'receipt_ids', at)
    for (const rid of receiptIds) {
      if (!(rid in receipts)) throw new BridgeContractError(`${at}: dangling receipt ${rid}`)
    }
    domains[id] = {
      kind: 'epistemic_domain', id,
      label: str(d, 'label', at),
      standpoint: str(d, 'standpoint', at),
      epi_status: oneOf(EPI_STATUSES, d['epi_status'], `${at}.epi_status`),
      receipt_ids: receiptIds,
    }
  }

  const nodes: BridgeIndex['nodes'] = {}
  for (const [id, n] of Object.entries(record(raw, 'nodes', 'index'))) {
    const at = `nodes.${id}`
    if (!isObj(n)) throw new BridgeContractError(`${at}: expected object`)
    const ds = strArr(n, 'domains', at)
    for (const did of ds) {
      if (!(did in domains)) throw new BridgeContractError(`${at}: dangling domain ${did}`)
    }
    nodes[id] = { label: str(n, 'label', at), rank: str(n, 'rank', at), domains: ds }
  }

  return {
    schema: BRIDGE_INDEX_SCHEMA,
    tree: { newick_sha256: str(tree, 'newick_sha256', 'index.tree'), node_count: nodeCount },
    nodes, domains, receipts,
  }
}

export function parseReceipt(raw: unknown, expectedId: string): ProofTransportReceipt {
  const at = `receipt(${expectedId})`
  if (!isObj(raw)) throw new BridgeContractError(`${at}: expected object`)
  if (raw['schema'] !== RECEIPT_SCHEMA) throw new BridgeContractError(`${at}.schema: expected ${RECEIPT_SCHEMA}`)
  const id = str(raw, 'id', at)
  // A receipt served under the wrong href is the "ReplayedArtifact" negative
  // case from epistemic-types, surfaced at the JSON layer.
  if (id !== expectedId) throw new BridgeContractError(`${at}: id mismatch (${id})`)
  const status = oneOf(PROOF_STATUSES, raw['status'], `${at}.status`)
  const agda = record(raw, 'agda', at)
  const safe = agda['safe']
  if (typeof safe !== 'boolean') throw new BridgeContractError(`${at}.agda.safe: expected boolean`)

  const out: ProofTransportReceipt = {
    schema: RECEIPT_SCHEMA, id,
    holder: str(raw, 'holder', at),
    artifact: str(raw, 'artifact', at),
    claim: str(raw, 'claim', at),
    meaning: str(raw, 'meaning', at),
    status,
    mode: oneOf(TRANSPORT_MODES, raw['mode'], `${at}.mode`),
    agda: {
      module: str(agda, 'module', `${at}.agda`),
      theorem: str(agda, 'theorem', `${at}.agda`),
      source: str(agda, 'source', `${at}.agda`),
      checked_commit: str(agda, 'checked_commit', `${at}.agda`),
      safe,
    },
  }
  const certified = status === 'Proof' || status === 'ProofUnder'
  if (raw['gap'] !== undefined) {
    if (certified) throw new BridgeContractError(`${at}: a ${status} cannot carry a gap`)
    out.gap = oneOf(GAPS, raw['gap'], `${at}.gap`)
  }
  // Mirrors opaqueNotCertifying: OpaqueReceipt mode can never be Proof.
  if (certified && out.mode === 'OpaqueReceipt') {
    throw new BridgeContractError(`${at}: OpaqueReceipt mode cannot be ${status}`)
  }
  if (status === 'ProofUnder') out.under = strArr(raw, 'under', at)
  if (isObj(raw['exacts'])) {
    const e = raw['exacts']
    out.exacts = {
      summary_ref: str(e, 'summary_ref', `${at}.exacts`),
      policy_fingerprint: str(e, 'policy_fingerprint', `${at}.exacts`),
    }
  }
  return out
}

// ── selection policy ─────────────────────────────────────────────────────────

/** Higher = shown first. Opinionated, and the only place the ordering lives. */
const STATUS_RANK: Record<ProofStatus, number> = {
  Proof: 5, ProofUnder: 4, Receipt: 3, Claimed: 2, Code: 1, Data: 0,
}

export interface Selection {
  nodeId: string
  label: string
  rank: string
  domains: EpistemicDomain[]
  /** The single receipt the panel renders, or null (nothing relevant). */
  chosen: ReceiptStub | null
  /** How many further relevant receipts exist — shown as a count, not rendered. */
  alternatives: number
}

/**
 * Click → selection, O(domains × receipts-per-domain) with no tree walk.
 * Returns null for ids not in the index (leaves, or a stale tree).
 *
 * Relevance: a receipt whose `artifacts` names this node beats a domain-wide
 * one (empty `artifacts`); receipts about OTHER nodes are not relevant at all.
 * Ties break on status strength, then id, so the result is deterministic.
 */
export function resolveSelection(index: BridgeIndex, nodeId: string): Selection | null {
  const node = index.nodes[nodeId]
  if (!node) return null
  const domains: EpistemicDomain[] = []
  const seen = new Set<string>()
  const scored: { stub: ReceiptStub; specific: boolean }[] = []
  for (const did of node.domains) {
    const d = index.domains[did]
    if (!d) continue
    domains.push(d)
    for (const rid of d.receipt_ids) {
      if (seen.has(rid)) continue
      seen.add(rid)
      const stub = index.receipts[rid]
      if (!stub) continue
      const specific = stub.artifacts.includes(nodeId)
      if (specific || stub.artifacts.length === 0) scored.push({ stub, specific })
    }
  }
  scored.sort((a, b) =>
    Number(b.specific) - Number(a.specific) ||
    STATUS_RANK[b.stub.status] - STATUS_RANK[a.stub.status] ||
    (a.stub.id < b.stub.id ? -1 : a.stub.id > b.stub.id ? 1 : 0))
  return {
    nodeId, label: node.label, rank: node.rank, domains,
    chosen: scored[0]?.stub ?? null,
    alternatives: Math.max(0, scored.length - 1),
  }
}
