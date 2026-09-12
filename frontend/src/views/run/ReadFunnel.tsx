// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useMemo, useState } from 'react'
import { useApi } from '../../hooks/useApi'
import { api } from '../../api/client'
import { PlotlyChart } from '../../components/PlotlyChart'

const pct = (a: number | null | undefined, b: number | null | undefined) =>
  a == null || !b ? '' : `${(100 * a / b).toFixed(1)}%`

export function ReadFunnel({ study, run, group }: { study: string; run: string; group?: string }) {
  const fetcher = useCallback(() => api.analysis.readFunnel(study, run, group), [study, run, group])
  const { data } = useApi(fetcher)
  const [open, setOpen] = useState(false)

  const figure = useMemo(() => {
    if (!data) return null
    const totals = data.stages.map(st =>
      data.samples.reduce((a, s) => a + (s.values[st.key] ?? 0), 0))
    return {
      data: [{
        type: 'funnel',
        y: data.stages.map(s => s.label),
        x: totals,
        textinfo: 'value+percent initial',
        marker: { color: '#4c9be8' },
      }],
      layout: { margin: { l: 140, r: 20, t: 10, b: 30 } },
    }
  }, [data])

  if (!data || !figure) return null
  const first = data.stages[0]?.key
  const last = data.stages[data.stages.length - 1]?.key

  return (
    <div className="card" style={{ marginTop: 16 }}>
      <div className="card-title">Reads through the pipeline</div>
      <PlotlyChart figure={figure} heightRatio={0.35} />
      {data.checks.length > 0 && (
        <ul style={{ listStyle: 'none', margin: '8px 0', fontSize: '.85rem' }}>
          {data.checks.map(c => (
            <li key={c.name} style={{ color: c.ok ? 'inherit' : 'var(--color-danger)' }}>
              <span aria-hidden="true">{c.ok ? '✓' : '✗'}</span> {c.name}: {c.detail}
            </li>
          ))}
        </ul>
      )}
      <button type="button" className="btn btn-sm" aria-expanded={open} onClick={() => setOpen(o => !o)}>
        {open ? 'Hide' : 'Show'} per-sample counts
      </button>
      {open && (
        <div style={{ overflowX: 'auto', marginTop: 8 }}>
          <table style={{ borderCollapse: 'collapse', fontSize: '.8rem' }}>
            <thead>
              <tr>
                <th style={{ textAlign: 'left', padding: '2px 8px' }}>Sample</th>
                {data.stages.map(s => <th key={s.key} style={{ textAlign: 'right', padding: '2px 8px' }}>{s.label}</th>)}
                <th style={{ textAlign: 'right', padding: '2px 8px' }}>Kept</th>
              </tr>
            </thead>
            <tbody>
              {data.samples.map(s => (
                <tr key={s.sample}>
                  <td style={{ padding: '2px 8px' }}>{s.sample}</td>
                  {data.stages.map(st => (
                    <td key={st.key} style={{ textAlign: 'right', padding: '2px 8px' }}>
                      {s.values[st.key] == null ? '' : s.values[st.key]!.toLocaleString()}
                    </td>
                  ))}
                  <td style={{ textAlign: 'right', padding: '2px 8px' }}>{pct(s.values[last], s.values[first])}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  )
}
