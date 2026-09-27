// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell; 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
//
// Builds the right-hand panel as a detached DocumentFragment. textContent only
// — receipt strings come from JSON and are never parsed as HTML.

import type { PanelState } from './controller'

export function renderPanel(s: PanelState, doc: Document): DocumentFragment {
  const frag = doc.createDocumentFragment()
  const el = <K extends keyof HTMLElementTagNameMap>(tag: K, text?: string, cls?: string) => {
    const e = doc.createElement(tag)
    if (text !== undefined) e.textContent = text
    if (cls) e.className = cls
    return e
  }
  if (s.kind === 'idle') { frag.append(el('p', 'Select an internal node.', 'pt-hint')); return frag }

  const { selection } = s
  frag.append(el('h3', `${selection.label} (${selection.rank})`))
  const domains = el('ul', undefined, 'pt-domains')
  for (const d of selection.domains) {
    const li = el('li', `${d.label} — ${d.epi_status} @ ${d.standpoint}`)
    li.dataset['epi'] = d.epi_status
    domains.append(li)
  }
  frag.append(domains)

  switch (s.kind) {
    case 'loading': frag.append(el('p', `Loading ${selection.chosen?.claim ?? ''}…`, 'pt-hint')); break
    case 'empty': frag.append(el('p', 'No ProofTransport receipt is relevant to this clade.', 'pt-hint')); break
    case 'error': frag.append(el('p', `Receipt unavailable: ${s.message}`, 'pt-error')); break
    case 'ready': {
      const r = s.receipt
      const card = el('article', undefined, 'pt-receipt')
      card.dataset['status'] = r.status
      const dl = el('dl')
      const row = (k: string, v: string) => dl.append(el('dt', k), el('dd', v))
      row('Status', r.gap ? `${r.status} (${r.gap})` : r.status)
      row('Claim', r.claim)
      row('Meaning', r.meaning)
      row('Holder → artefact', `${r.holder} → ${r.artifact}`)
      row('Mode', r.mode)
      if (r.under?.length) row('Under', r.under.join('; '))
      row('Agda', `${r.agda.module}.${r.agda.theorem}${r.agda.safe ? ' (--safe)' : ''}`)
      row('Checked at', r.agda.checked_commit.slice(0, 12))
      if (r.exacts) row('Exacts', `${r.exacts.summary_ref} [${r.exacts.policy_fingerprint.slice(0, 12)}]`)
      card.append(dl, el('p', 'Displayed, not re-verified in the browser.', 'pt-hint'))
      frag.append(card)
      if (selection.alternatives > 0) {
        frag.append(el('p', `${selection.alternatives} further relevant receipt(s) not shown.`, 'pt-hint'))
      }
    }
  }
  return frag
}
