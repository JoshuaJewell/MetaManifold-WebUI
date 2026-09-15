// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useState } from 'react'
import { api } from '../api/client'
import { AnalysisChart } from './AnalysisChart'
import { useAlphaMetricSelection, AlphaMetricToggles, ALPHA_METRICS, extractAlphaPanel } from './alphaMetrics'
import { useToast } from './Toast'
import { errorMessage } from '../api/errorMessage'
import type { AnalysisOption } from './annotationShared'
import type { ColFilter, ComparisonRunSpec, PermanovaResult } from '../api/types'

/** Alpha comparison, NMDS and PERMANOVA across the runs/sub-groups in `runs`. */
export function DiversityPanel({ study, runs, option, aggregate }: {
  study: string
  runs: ComparisonRunSpec[]
  option: AnalysisOption | null
  aggregate: boolean
}) {
  const toast = useToast()
  const [alphaFig, setAlphaFig]   = useState<unknown>(null)
  const [nmdsFig, setNmdsFig]     = useState<unknown>(null)
  const [permanova, setPermanova] = useState<PermanovaResult | null>(null)
  const [loading, setLoading]     = useState<string | null>(null)
  const [rAvailable, setRAvailable] = useState<boolean | null>(null)

  useEffect(() => {
    api.analysis.capabilities().then(c => setRAvailable(c.r_available)).catch(() => setRAvailable(false))
  }, [])

  const body = option && runs.length >= 2
    ? {
        table: option.table,
        runs,
        colFilters: {} as Record<string, ColFilter>,
        aggregate,
      }
    : null

  const run = async (type: 'alpha' | 'nmds' | 'permanova') => {
    if (!body) return
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
      const label = type === 'alpha' ? 'Alpha comparison' : type === 'nmds' ? 'NMDS' : 'PERMANOVA'
      toast.error(`${label} failed: ${errorMessage(err)}`)
    } finally {
      setLoading(null)
    }
  }

  const { metrics, toggle } = useAlphaMetricSelection()

  if (runs.length < 2) {
    return <p className="empty-state">Select at least two runs or sub-groups to compare diversity.</p>
  }

  return (
    <div>
      <div style={{ display: 'flex', gap: 8, marginBottom: 12, flexWrap: 'wrap', alignItems: 'center' }}>
        <button className="btn" onClick={() => run('alpha')} disabled={loading !== null || !body}>
          {loading === 'alpha' ? 'Computing…' : 'Alpha Comparison'}
        </button>
        {rAvailable && (
          <button className="btn" onClick={() => run('nmds')} disabled={loading !== null || !body}>
            {loading === 'nmds' ? 'Computing…' : 'NMDS'}
          </button>
        )}
        {rAvailable && (
          <button className="btn" onClick={() => run('permanova')} disabled={loading !== null || !body}>
            {loading === 'permanova' ? 'Computing…' : 'PERMANOVA'}
          </button>
        )}
        {rAvailable === false && (
          <span style={{ fontSize: '.82rem', color: 'var(--color-muted-fg)' }}>
            NMDS and PERMANOVA need R with vegan.
          </span>
        )}
      </div>

      {alphaFig != null && (
        <>
          <AlphaMetricToggles metrics={metrics} toggle={toggle} />
          {ALPHA_METRICS.filter(m => metrics.has(m)).map(m => (
            <AnalysisChart key={m} study={study} figure={extractAlphaPanel(alphaFig, m)} heightRatio={0.32} />
          ))}
        </>
      )}
      {nmdsFig != null && (
        <AnalysisChart study={study} figure={nmdsFig} heightRatio={0.48} />
      )}
      {permanova && (
        <div className="card" style={{ fontFamily: 'monospace', fontSize: '.82rem', whiteSpace: 'pre-wrap' }}>
          <div className="card-title">PERMANOVA Results</div>
          <pre style={{ margin: 0 }}>{permanova.text}</pre>
          {(permanova.terms?.length
            ? permanova.terms
            : permanova.p_value != null
              ? [{ term: 'model', r2: permanova.r2 ?? NaN, f_statistic: permanova.f_statistic, p_value: permanova.p_value }]
              : []
          ).map(t => (
            <p key={t.term} style={{ marginTop: 8, fontFamily: 'inherit' }}>
              {t.term}: R&sup2; = {t.r2.toFixed(3)} · F = {t.f_statistic?.toFixed(2) ?? 'n/a'} · p = {t.p_value?.toFixed(4) ?? 'n/a'}
            </p>
          ))}
          {(permanova.terms?.length ?? 0) > 1 && (
            <p style={{ marginTop: 4, fontFamily: 'inherit' }}>Terms are tested in order, each after the ones above it.</p>
          )}
          {permanova.blocked && (
            <p style={{ marginTop: 4, fontFamily: 'inherit' }}>Permutations restricted within individuals.</p>
          )}
          {permanova.dispersion && (
            <>
              <div className="card-title" style={{ marginTop: 16 }}>PERMDISP (homogeneity of dispersion)</div>
              {'error' in permanova.dispersion ? (
                <p style={{ fontFamily: 'inherit' }}>{permanova.dispersion.error}</p>
              ) : (
                <>
                  <pre style={{ margin: 0 }}>{permanova.dispersion.text}</pre>
                  <p style={{ marginTop: 8, fontFamily: 'inherit' }}>
                    F({permanova.dispersion.df[0]}, {permanova.dispersion.df[1]}) = {permanova.dispersion.f_statistic.toFixed(2)} · p = {permanova.dispersion.p_value.toFixed(4)}
                  </p>
                  {permanova.dispersion.groups.map(g => (
                    <p key={g.group} style={{ marginTop: 4, fontFamily: 'inherit' }}>
                      {g.group}: distance to centroid, mean {g.mean.toFixed(3)} · median {g.median.toFixed(3)}
                    </p>
                  ))}
                </>
              )}
            </>
          )}
        </div>
      )}
    </div>
  )
}
