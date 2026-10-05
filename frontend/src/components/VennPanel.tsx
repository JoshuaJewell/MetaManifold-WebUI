// frontend/src/components/VennPanel.tsx
// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).

import { useState, useEffect, useMemo, useRef } from 'react'
import { VennDiagram, UpSetJS, asSets } from '@upsetjs/react'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import { useToast } from './Toast'
import type { ComparisonRunSpec, VennResult } from '../api/types'
import { sharedRanks } from './sharedRanks'
import type { AnalysisOption } from './annotationShared'

type Mode = 'euler' | 'upset'

/** Taxon overlap between the selected runs at one shared rank, as an Euler or UpSet diagram. */
export function VennPanel({ study, runs, option }: {
  study: string
  runs: ComparisonRunSpec[]
  option: AnalysisOption | null
}) {
  const toast = useToast()
  const [ranks, setRanks]   = useState<string[]>([])
  const [rank, setRank]     = useState<string>('')
  const [mode, setMode]     = useState<Mode>('euler')
  const [result, setResult] = useState<VennResult | null>(null)
  const [loading, setLoading] = useState(false)
  const [ranksReady, setRanksReady] = useState(false)
  const [ranksError, setRanksError] = useState<string | null>(null)
  // Only the response to the latest request, for the current inputs, is shown.
  const request = useRef(0)
  const box = useRef<HTMLDivElement>(null)
  const [width, setWidth] = useState(580)

  useEffect(() => {
    const el = box.current
    if (!el) return
    const ro = new ResizeObserver(entries => setWidth(Math.max(240, Math.min(580, entries[0].contentRect.width))))
    ro.observe(el)
    return () => ro.disconnect()
  }, [])

  // Discover the intersection of taxonomy ranks available across all selected runs.
  useEffect(() => {
    setRanksError(null)
    if (!option || runs.length === 0) { setRanks([]); setRank(''); setRanksReady(false); return }
    let cancelled = false
    setRanksReady(false)
    sharedRanks(runs, r => api.analysis.ranks(study, r.run, {
      table: option.table,
      group: r.group,
      ...(r.source ? { source: r.source } : {}),
    })).then(intersection => {
      if (cancelled) return
      setRanks(intersection)
      setRanksReady(true)
      // Default to the deepest (last) shared rank.
      setRank(current =>
        intersection.includes(current) ? current : (intersection[intersection.length - 1] ?? '')
      )
    }).catch(err => {
      if (cancelled) return
      setRanks([])
      setRank('')
      setRanksError(errorMessage(err))
    })
    return () => { cancelled = true }
  }, [study, runs, option])

  // A diagram answers the inputs it was computed from; new inputs discard it,
  // along with any request still in flight for the old ones.
  useEffect(() => {
    request.current++
    setResult(null)
    setLoading(false)
  }, [study, runs, option?.table, rank])

  /** Compute the overlap for the current inputs and show it, unless the inputs changed meanwhile. */
  const runVenn = async () => {
    if (!option || !rank) return
    const req = ++request.current
    setLoading(true)
    try {
      const res = await api.analysis.venn(study, {
        runs,
        table: option.table,
        rank,
      })
      if (req === request.current) setResult(res)
    } catch (err) {
      if (req !== request.current) return
      setResult(null)
      toast.error(`Taxon overlap failed: ${errorMessage(err)}`)
    } finally {
      if (req === request.current) setLoading(false)
    }
  }

  const upsetjsSets = useMemo(
    () => result ? asSets(result.sets.map(s => ({ name: s.name, elems: s.taxa }))) : null,
    [result],
  )

  return (
    <div style={{ marginTop: 16 }}>
      <div style={{ display: 'flex', gap: 8, alignItems: 'center', flexWrap: 'wrap', marginBottom: 12 }}>
        <button
          className="btn"
          onClick={runVenn}
          disabled={loading || !option || !rank}
        >
          {loading ? 'Computing…' : 'Taxon Overlap'}
        </button>

        {ranks.length > 0 && (
          <select
            aria-label="Rank"
            value={rank}
            onChange={e => setRank(e.target.value)}
            style={{ font: 'inherit', padding: '1px 4px', verticalAlign: 'baseline' }}
          >
            {ranks.map(r => <option key={r} value={r}>{r}</option>)}
          </select>
        )}
        {ranksReady && ranks.length === 0 && (
          <span style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)' }}>No taxonomy rank shared by the selected runs.</span>
        )}
        {ranksError && (
          <span className="error-msg" style={{ margin: 0 }}>Could not read the taxonomy ranks: {ranksError}</span>
        )}

        <div style={{ display: 'flex', gap: 4 }} role="group" aria-label="Diagram type">
          <button className={`btn btn-sm ${mode === 'euler' ? 'btn-primary' : ''}`} aria-pressed={mode === 'euler'}
            onClick={() => setMode('euler')}>Euler</button>
          <button className={`btn btn-sm ${mode === 'upset' ? 'btn-primary' : ''}`} aria-pressed={mode === 'upset'}
            onClick={() => setMode('upset')}>UpSet</button>
        </div>
      </div>

      <div ref={box} style={{ maxWidth: 580 }}>
        {upsetjsSets && (mode === 'euler'
          ? <VennDiagram sets={upsetjsSets} width={width} height={Math.round(width * 340 / 580)} />
          : <UpSetJS sets={upsetjsSets} width={width} height={Math.round(width * 340 / 580)} />)}
      </div>
    </div>
  )
}
