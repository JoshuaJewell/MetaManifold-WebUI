// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
// type-theory section: dendrogram → bridge index → ProofTransport receipt.
// DOM-free: exercises the parser, the selection policy and the controller's
// abort / stale-response / cache behaviour against the checked-in examples.
import { test, expect } from 'bun:test'
import { readFileSync } from 'node:fs'
import { join } from 'node:path'
import {
  BridgeContractError, parseBridgeIndex, parseReceipt, resolveSelection,
} from '../../src/type-theory/dendrogram-receipts/bridgeIndex'
import { createReceiptController, Lru, type PanelState } from '../../src/type-theory/dendrogram-receipts/controller'

const EX = join(import.meta.dir, '../../../type-theory/dendrogram-receipts/examples')
const load = (p: string): unknown => JSON.parse(readFileSync(join(EX, p), 'utf8'))
const index = parseBridgeIndex(load('bridge-index.example.json'))

test('every example receipt parses under its own id', () => {
  for (const [id, stub] of Object.entries(index.receipts)) {
    expect(parseReceipt(load(stub.href), id).id).toBe(id)
  }
})

test('node-specific receipt beats a stronger-looking domain-wide one', () => {
  const s = resolveSelection(index, 'n2')!
  expect(s.domains.map((d) => d.id)).toEqual(['exacts'])
  expect(s.chosen?.id).toBe('r-checked-sum-n2')
  expect(s.alternatives).toBe(1)
})

test('without a specific receipt, status strength decides; other nodes\' receipts are excluded', () => {
  const s = resolveSelection(index, 'n1')!
  expect(s.chosen?.id).toBe('r-checked-sum-any') // ProofUnder > Receipt; n2-specific excluded
  expect(s.alternatives).toBe(1)
})

test('leaves / unknown ids resolve to null (click is a no-op)', () => {
  expect(resolveSelection(index, 'leaf-17')).toBeNull()
})

test('contract rejections', () => {
  const bad = load('bridge-index.example.json') as { receipts: Record<string, { href: string }> }
  bad.receipts['r-checked-sum-n2']!.href = 'https://evil.example/r.json'
  expect(() => parseBridgeIndex(bad)).toThrow(BridgeContractError)

  const r = load('receipts/r-epi-not-factive.json') as Record<string, unknown>
  expect(() => parseReceipt(r, 'someone-else')).toThrow(/id mismatch/)
  expect(() => parseReceipt({ ...r, status: 'Proof' }, 'r-epi-not-factive')).toThrow(/gap/)
  const { gap: _gap, ...noGap } = r
  expect(() => parseReceipt({ ...noGap, status: 'Proof' }, 'r-epi-not-factive')).toThrow(/OpaqueReceipt/)
})

test('LRU evicts least-recently-used', () => {
  const l = new Lru<string, number>(2)
  l.set('a', 1); l.set('b', 2); l.get('a'); l.set('c', 3)
  expect(l.get('b')).toBeUndefined()
  expect(l.get('a')).toBe(1)
})

test('controller: rapid clicks abort the old fetch and never render a stale receipt', async () => {
  const states: PanelState[] = []
  const aborted: string[] = []
  const gates = new Map<string, () => void>()
  const ctl = createReceiptController({
    index,
    onState: (s) => states.push(s),
    fetchJson: (href, signal) => new Promise((res) => {
      signal.addEventListener('abort', () => aborted.push(href))
      gates.set(href, () => res(load(href)))
    }),
  })
  const p1 = ctl.select('n2')
  const p2 = ctl.select('n1')
  gates.get('receipts/r-checked-sum-n2.json')!() // late response for the first click
  gates.get('receipts/r-checked-sum-any.json')!()
  await Promise.all([p1, p2])

  expect(aborted).toEqual(['receipts/r-checked-sum-n2.json'])
  const ready = states.filter((s) => s.kind === 'ready')
  expect(ready.length).toBe(1)
  expect(ready[0]!.kind === 'ready' && ready[0]!.receipt.id).toBe('r-checked-sum-any')
  expect(ctl.cacheSize).toBe(2) // the stale one is still cached for next time

  // Returning to n2 is served from cache: no loading state, no fetch.
  const before = states.length
  await ctl.select('n2')
  expect(states.slice(before).map((s) => s.kind)).toEqual(['ready'])
})

test('controller: re-clicking the same node does nothing', async () => {
  let fetches = 0
  const ctl = createReceiptController({
    index, onState: () => {},
    fetchJson: async (href) => { fetches++; return load(href) },
  })
  await ctl.select('n2'); await ctl.select('n2')
  expect(fetches).toBe(1)
})
