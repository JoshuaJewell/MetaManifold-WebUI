// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell; 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
//
// Click → bridge index → one receipt → panel, without choking the DOM.
//
// The five rules that keep it cheap (see DESIGN.md §"Not choking the DOM"):
//   1. ONE delegated listener on the tree root, never one per node.
//   2. The click handler does no layout reads and no fetch-waiting: it resolves
//      the selection synchronously from the in-memory index (O(1) lookup).
//   3. Every new selection aborts the previous fetch; a generation counter
//      discards any response that still races in.
//   4. Full receipts are fetched lazily and kept in a small LRU.
//   5. The panel is written at most once per animation frame, as a single
//      replaceChildren() of a detached fragment built with textContent only.

import { parseReceipt, resolveSelection, type Selection } from './bridgeIndex'
import type { BridgeIndex, DendrogramHost, ProofTransportReceipt } from './types'
import { renderPanel } from './renderPanel'

export type PanelState =
  | { kind: 'idle' }
  | { kind: 'loading'; selection: Selection }
  | { kind: 'empty'; selection: Selection }
  | { kind: 'ready'; selection: Selection; receipt: ProofTransportReceipt }
  | { kind: 'error'; selection: Selection; message: string }

export class Lru<K, V> {
  private readonly m = new Map<K, V>()
  constructor(private readonly cap: number) {}
  get(k: K): V | undefined {
    const v = this.m.get(k)
    if (v !== undefined) { this.m.delete(k); this.m.set(k, v) }
    return v
  }
  set(k: K, v: V): void {
    this.m.delete(k)
    this.m.set(k, v)
    if (this.m.size > this.cap) this.m.delete(this.m.keys().next().value as K)
  }
  get size(): number { return this.m.size }
}

export interface ControllerDeps {
  index: BridgeIndex
  fetchJson: DendrogramHost['fetchJson']
  /** Receives every state; the DOM sink coalesces, tests just record. */
  onState: (s: PanelState) => void
  cacheSize?: number
}

/** DOM-free state machine. `select` is safe to call at pointer-event rate. */
export function createReceiptController(deps: ControllerDeps) {
  const cache = new Lru<string, ProofTransportReceipt>(deps.cacheSize ?? 32)
  let generation = 0
  let inflight: AbortController | null = null
  let current: string | null = null

  async function select(nodeId: string): Promise<void> {
    if (nodeId === current) return // re-click on the same node: no work
    const selection = resolveSelection(deps.index, nodeId)
    if (!selection) return // leaf / unknown id: keep the current panel
    current = nodeId
    const gen = ++generation
    inflight?.abort()
    inflight = null

    const stub = selection.chosen
    if (!stub) { deps.onState({ kind: 'empty', selection }); return }
    const hit = cache.get(stub.id)
    if (hit) { deps.onState({ kind: 'ready', selection, receipt: hit }); return }

    deps.onState({ kind: 'loading', selection })
    const ac = new AbortController()
    inflight = ac
    try {
      const receipt = parseReceipt(await deps.fetchJson(stub.href, ac.signal), stub.id)
      cache.set(stub.id, receipt) // cache even if stale: the next click may want it
      if (gen !== generation) return
      deps.onState({ kind: 'ready', selection, receipt })
    } catch (e) {
      if (gen !== generation || ac.signal.aborted) return
      deps.onState({ kind: 'error', selection, message: e instanceof Error ? e.message : String(e) })
    } finally {
      if (inflight === ac) inflight = null
    }
  }

  function dispose(): void {
    generation++
    inflight?.abort()
    inflight = null
  }

  return { select, dispose, get cacheSize() { return cache.size } }
}

/**
 * Wire the controller to a real DOM. The renderer must mark internal nodes as
 *   <g data-node-id="n42" data-node-kind="internal" tabindex="0">…</g>
 * (Protoctist's NHX `ND=` id is the natural source of data-node-id.)
 * Returns a disposer; call it when the tree is re-rendered or unmounted.
 */
export function attachDendrogram(host: DendrogramHost, index: BridgeIndex): () => void {
  let pending: PanelState | null = null
  let frame = 0
  let selected: Element | null = null
  const flush = () => {
    frame = 0
    if (pending) { host.panel.replaceChildren(renderPanel(pending, host.panel.ownerDocument)); pending = null }
  }
  const ctl = createReceiptController({
    index,
    fetchJson: host.fetchJson,
    onState: (s) => { pending = s; if (!frame) frame = requestAnimationFrame(flush) },
  })

  const nodeFrom = (t: EventTarget | null): Element | null =>
    t instanceof Element ? t.closest('[data-node-kind="internal"][data-node-id]') : null

  const activate = (el: Element) => {
    // Selection highlight is one attribute flip on two elements, not a re-render.
    selected?.removeAttribute('aria-selected')
    el.setAttribute('aria-selected', 'true')
    selected = el
    void ctl.select(el.getAttribute('data-node-id') ?? '')
  }
  const onClick = (ev: Event) => { const el = nodeFrom(ev.target); if (el) activate(el) }
  const onKey = (ev: Event) => {
    const k = (ev as KeyboardEvent).key
    if (k !== 'Enter' && k !== ' ') return
    const el = nodeFrom(ev.target)
    if (el) { ev.preventDefault(); activate(el) }
  }

  host.panel.setAttribute('aria-live', 'polite')
  host.treeRoot.addEventListener('click', onClick)
  host.treeRoot.addEventListener('keydown', onKey)
  return () => {
    host.treeRoot.removeEventListener('click', onClick)
    host.treeRoot.removeEventListener('keydown', onKey)
    if (frame) cancelAnimationFrame(frame)
    ctl.dispose()
  }
}
