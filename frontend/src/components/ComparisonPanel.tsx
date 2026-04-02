import { useEffect, useState } from 'react'
import { api } from '../api/client'
import { PlotlyChart } from './PlotlyChart'
import { useToast } from './Toast'
import type { ComparisonRunSpec, ColFilter, PermanovaResult } from '../api/types'

export function ComparisonPanel({ study, runs }: {
  study: string
  runs: ComparisonRunSpec[]
}) {
  const toast = useToast()
  const [alphaFig, setAlphaFig]       = useState<unknown>(null)
  const [nmdsFig, setNmdsFig]         = useState<unknown>(null)
  const [permanova, setPermanova]     = useState<PermanovaResult | null>(null)
  const [loading, setLoading]         = useState<string | null>(null)
  const [tables, setTables]           = useState<string[]>([])
  const [table, setTable]             = useState('merged')
  const [rAvailable, setRAvailable]   = useState<boolean | null>(null)

  useEffect(() => {
    api.analysis.capabilities().then(c => setRAvailable(c.r_available)).catch(() => setRAvailable(false))
  }, [])

  const isSubgroupMode = runs.length > 0 && runs.every(r => r.prefix)

  // Discover tables available across the selected runs (deduplicate by run+group)
  useEffect(() => {
    let cancelled = false
    const seen = new Set<string>()
    const unique = runs.filter(r => {
      const key = `${r.run}|${r.group ?? ''}`
      if (seen.has(key)) return false
      seen.add(key)
      return true
    })
    Promise.all(
      unique.map(r => api.results.runTables(study, r.run, r.group).catch(() => []))
    ).then(results => {
      if (cancelled) return
      // Union of all table IDs
      const ids = new Set<string>()
      for (const metas of results) for (const m of metas) ids.add(m.id)
      const sorted = [...ids].sort((a, b) => a === 'merged' ? -1 : b === 'merged' ? 1 : a.localeCompare(b))
      setTables(sorted)
      if (sorted.length > 0 && !sorted.includes(table)) setTable(sorted[0])
    })
    return () => { cancelled = true }
  }, [study, runs, table])

  const body = { table, runs, colFilters: {} as Record<string, ColFilter> }

  const run = async (type: 'alpha' | 'nmds' | 'permanova') => {
    setLoading(type)
    try {
      if (type === 'alpha') {
        setAlphaFig(await api.analysis.compareAlpha(study, body))
      } else if (type === 'nmds') {
        setNmdsFig(await api.analysis.nmds(study, body))
      } else {
        setPermanova(await api.analysis.permanova(study, body))
      }
    } catch (err) {
      toast.error(`${type} failed: ${err instanceof Error ? err.message : err}`)
    } finally {
      setLoading(null)
    }
  }

  if (runs.length < 2) return null

  return (
    <div style={{ marginTop: 24 }}>
      <h2 style={{ fontSize: '1rem', marginBottom: 8 }}>
        {isSubgroupMode ? 'Cross-Group Comparison' : 'Cross-Run Comparison'}
      </h2>
      <p style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)', marginBottom: 12 }}>
        Comparing {runs.length} {isSubgroupMode ? 'groups' : 'runs'} on the{' '}
        {tables.length > 1 ? (
          <select
            value={table}
            onChange={e => setTable(e.target.value)}
            style={{ font: 'inherit', padding: '1px 4px', verticalAlign: 'baseline' }}
          >
            {tables.map(t => <option key={t} value={t}>{t}</option>)}
          </select>
        ) : (
          <code>{table}</code>
        )}
        {' '}table
      </p>
      <div style={{ display: 'flex', gap: 8, marginBottom: 12 }}>
        <button className="btn" onClick={() => run('alpha')} disabled={loading !== null}>
          {loading === 'alpha' ? 'Computing...' : 'Alpha Comparison'}
        </button>
        {rAvailable && (
          <button className="btn" onClick={() => run('nmds')} disabled={loading !== null}>
            {loading === 'nmds' ? 'Computing...' : 'NMDS'}
          </button>
        )}
        {rAvailable && (
          <button className="btn" onClick={() => run('permanova')} disabled={loading !== null}>
            {loading === 'permanova' ? 'Computing...' : 'PERMANOVA'}
          </button>
        )}
      </div>

      {alphaFig != null && <PlotlyChart figure={alphaFig} />}
      {nmdsFig != null && <PlotlyChart figure={nmdsFig} />}
      {permanova && (
        <div className="card" style={{ fontFamily: 'monospace', fontSize: '.82rem', whiteSpace: 'pre-wrap' }}>
          <div className="card-title">PERMANOVA Results</div>
          <pre style={{ margin: 0 }}>{permanova.text}</pre>
          {permanova.p_value != null && (
            <p style={{ marginTop: 8, fontFamily: 'inherit' }}>
              R&sup2; = {permanova.r2?.toFixed(3)} | F = {permanova.f_statistic?.toFixed(2)} | p = {permanova.p_value?.toFixed(4)}
            </p>
          )}
        </div>
      )}
    </div>
  )
}
