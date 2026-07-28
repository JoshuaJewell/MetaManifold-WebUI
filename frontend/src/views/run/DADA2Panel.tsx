// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState, useEffect } from 'react'
import { api, apiUrl } from '../../api/client'
import { StageConfig } from '../../components/PipelineStages'
import type { ConfigMap, RunStages } from '../../api/types'

interface DADA2Output {
  figures: { name: string; label: string; url: string }[]
  has_stats: boolean
  logs: { name: string; url: string }[]
  config: Record<string, unknown>
}

type DADA2SubTab = 'quality' | 'denoising' | 'merging' | 'taxonomy'

const DADA2_SUBTABS: { key: DADA2SubTab; label: string; runTarget: string }[] = [
  { key: 'quality',   label: 'Quality & Filtering',  runTarget: 'filter_trim' },
  { key: 'denoising', label: 'Denoising',            runTarget: 'learn_errors' },
  { key: 'merging',   label: 'Merging & Lengths',    runTarget: 'denoise' },
  { key: 'taxonomy',  label: 'Taxonomy',             runTarget: 'assign_taxonomy' },
]

const SUBTAB_FIGURES: Record<DADA2SubTab, string[]> = {
  quality:   ['quality_unfiltered', 'quality_filtered'],
  denoising: ['error_rates'],
  merging:   ['length_distribution', 'length_distribution_filtered'],
  taxonomy:  [],
}

const SUBTAB_CONFIG: Record<DADA2SubTab, { prefix: string; label: string }[]> = {
  quality:   [{ prefix: 'dada2.filter_trim.', label: 'Filter & Trim' }],
  denoising: [{ prefix: 'dada2.dada.', label: 'Denoising' }],
  merging:   [
    { prefix: 'dada2.merge.', label: 'Merge Pairs' },
    { prefix: 'dada2.asv.', label: 'ASV & Chimera Filtering' },
  ],
  taxonomy:  [{ prefix: 'dada2.taxonomy.', label: 'Taxonomy Assignment' }],
}

const SUBTAB_HINTS: Record<DADA2SubTab, string> = {
  quality:   'Inspect pre-filter quality to choose truncation lengths, then verify post-filter output.',
  denoising: 'Error rate convergence plots. Re-run if denoising parameters change.',
  merging:   'Check length distribution to set band-size min/max, then verify filtering removed off-target amplicons.',
  taxonomy:  'Review pipeline stats for unexpected read loss before running taxonomy assignment.',
}

const SUBTAB_STALE_SECTIONS: Record<DADA2SubTab, string[]> = {
  quality:   ['dada2.file_patterns', 'dada2.filter_trim'],
  denoising: ['dada2.dada'],
  merging:   ['dada2.merge', 'dada2.asv'],
  taxonomy:  ['dada2.taxonomy', 'dada2.output'],
}

function isSubTabStale(tab: DADA2SubTab, staleKeys: string[]): boolean {
  if (staleKeys.length === 0) return false
  const sections = SUBTAB_STALE_SECTIONS[tab]
  return sections.some(sec => staleKeys.some(k => k.startsWith(sec)))
}

export function StaleKeysBadge({ staleKeys, configMap }: { staleKeys: string[]; configMap?: ConfigMap | null }) {
  if (!Array.isArray(staleKeys) || staleKeys.length === 0) return null

  const lines = staleKeys.map(k => {
    const entry = configMap?.[k]
    if (entry) {
      const v = typeof entry.value === 'object' ? JSON.stringify(entry.value) : String(entry.value)
      return `${k} = ${v} (${entry.source})`
    }
    const matchingKeys = configMap
      ? Object.entries(configMap).filter(([ck]) => ck.startsWith(k + '.'))
      : []
    if (matchingKeys.length > 0) {
      return matchingKeys.map(([ck, { value, source }]) => {
        const v = typeof value === 'object' ? JSON.stringify(value) : String(value)
        return `${ck} = ${v} (${source})`
      }).join('\n')
    }
    return k
  })

  return (
    <span
      style={{ fontSize: '.75rem', color: '#f59e0b', fontWeight: 500, cursor: 'default', position: 'relative' }}
      title={lines.join('\n')}
    >
      Config changed
    </span>
  )
}



export function DADA2Panel({ study, run, group, dada2Data, configMap, onConfigChanged, stages, onRunStage, cacheKey }: {
  study: string; run: string; group?: string; dada2Data: DADA2Output | null
  configMap: ConfigMap | null; onConfigChanged: () => void
  stages: RunStages | null; onRunStage: (stage: string) => void
  cacheKey?: string | null
}) {
  const [showStats, setShowStats] = useState(false)
  const [statsData, setStatsData] = useState<{ columns: string[]; rows: Record<string, unknown>[] } | null>(null)
  const [expandedLog, setExpandedLog] = useState<string | null>(null)
  const [logContent, setLogContent] = useState<Record<string, string>>({})

  useEffect(() => {
    setStatsData(null)
    setExpandedLog(null)
    setLogContent({})
  }, [cacheKey])

  // Auto-select earliest stale sub-tab, or default to 'quality'
  const denoiseStaleKeys = stages?.dada2_denoise?.stale_keys ?? []
  const classifyStaleKeys = stages?.dada2_classify?.stale_keys ?? []
  const allStaleKeys = [...denoiseStaleKeys, ...classifyStaleKeys]

  const initialTab = (): DADA2SubTab => {
    for (const { key } of DADA2_SUBTABS) {
      if (isSubTabStale(key, allStaleKeys)) return key
    }
    return 'quality'
  }
  const [subTab, setSubTab] = useState<DADA2SubTab>(initialTab)

  useEffect(() => {
    if (showStats && !statsData) {
      api.results.dada2Stats(study, run, group).then(setStatsData).catch(() => {})
    }
  }, [showStats, statsData, study, run, group])

  const fetchLog = async (url: string, name: string) => {
    if (logContent[name]) { setExpandedLog(expandedLog === name ? null : name); return }
    try {
      const res = await fetch(apiUrl(url))
      if (!res.ok) throw new Error(res.statusText)
      const text = await res.text()
      setLogContent(prev => ({ ...prev, [name]: text }))
      setExpandedLog(name)
    } catch {
      setLogContent(prev => ({ ...prev, [name]: 'Failed to load log.' }))
      setExpandedLog(name)
    }
  }

  const denoiseStatus  = stages?.dada2_denoise?.status ?? ''
  const denoiseRunning = denoiseStatus === 'running'
  const classifyStatus = stages?.dada2_classify?.status ?? ''
  const classifyRunning = classifyStatus === 'running'
  const isRunning = denoiseRunning || classifyRunning

  const tabDef = DADA2_SUBTABS.find(t => t.key === subTab)!
  const figNames = SUBTAB_FIGURES[subTab]
  const figures = figNames.map(name => dada2Data?.figures.find(f => f.name === name) ?? null)
  const configSections = SUBTAB_CONFIG[subTab]
  const hint = SUBTAB_HINTS[subTab]
  const tabStale = isSubTabStale(subTab, allStaleKeys)

  return (
    <>
      {/* Sub-tabs */}
      <div className="tabs" style={{ marginBottom: 0 }}>
        {DADA2_SUBTABS.map(t => {
          const stale = isSubTabStale(t.key, allStaleKeys)
          return (
            <button key={t.key} className={`tab ${subTab === t.key ? 'active' : ''}`}
              onClick={() => setSubTab(t.key)}>
              {t.label}
              {stale && <span style={{ display: 'inline-block', width: 7, height: 7, borderRadius: '50%', background: '#f59e0b', marginLeft: 6, verticalAlign: 'middle' }} title="Config changed" />}
            </button>
          )
        })}
      </div>

      {/* Panel content */}
      <div className="card" style={{ borderTopLeftRadius: 0, borderTopRightRadius: 0 }}>
        {/* Figures: side-by-side for pairs, single for solo */}
        {figNames.length > 0 && (
          <div style={{ display: 'flex', flexWrap: 'wrap', gap: 12, marginBottom: 12 }}>
            {figures.map((fig, i) => (
              <div key={figNames[i]} style={{ flex: 1, minWidth: 320 }}>
                <div style={{ fontSize: '.78rem', fontWeight: 600, color: 'var(--color-muted-fg)', marginBottom: 4 }}>
                  {fig?.label ?? figNames[i].replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase())}
                </div>
                {fig ? (
                  <iframe
                    src={apiUrl(fig.url) + (cacheKey ? `?v=${encodeURIComponent(cacheKey)}` : '')}
                    style={{ width: '100%', height: figures.length > 1 ? 'calc(50vh - 120px)' : 'calc(100vh - 380px)', minHeight: 280, border: '1px solid var(--color-border)', borderRadius: 4 }}
                    title={fig.label}
                  />
                ) : (
                  <div style={{ height: 200, display: 'flex', alignItems: 'center', justifyContent: 'center', background: 'var(--color-surface)', border: '1px solid var(--color-border)', borderRadius: 4, color: 'var(--color-muted-fg)', fontSize: '.85rem' }}>
                    Not yet generated
                  </div>
                )}
              </div>
            ))}
          </div>
        )}

        {/* Taxonomy: pipeline stats instead of figures */}
        {subTab === 'taxonomy' && (
          <>
            <button type="button" aria-expanded={showStats} disabled={!dada2Data?.has_stats}
              style={{ fontSize: '.85rem', fontWeight: 600, marginBottom: 8, border: 'none', background: 'none', padding: 0,
                       font: 'inherit', color: 'inherit', cursor: dada2Data?.has_stats ? 'pointer' : 'default' }}
              onClick={() => setShowStats(!showStats)}>
              Pipeline Stats
              {dada2Data?.has_stats && (
                <span style={{ fontSize: '.78rem', color: 'var(--color-muted-fg)', marginLeft: 8, fontWeight: 400 }}>{showStats ? 'Hide' : 'Show'}</span>
              )}
            </button>
            {dada2Data?.has_stats && showStats && (
              statsData ? (
                <div style={{ overflowX: 'auto', marginBottom: 12 }}>
                  <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '.82rem' }}>
                    <thead>
                      <tr>
                        {statsData.columns.map(c => (
                          <th key={c} style={{ padding: '6px 10px', borderBottom: '2px solid var(--color-border)', textAlign: 'left', whiteSpace: 'nowrap' }}>{c}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {statsData.rows.map((row, i) => (
                        <tr key={i}>
                          {statsData.columns.map(c => (
                            <td key={c} style={{ padding: '4px 10px', borderBottom: '1px solid var(--color-border)', whiteSpace: 'nowrap' }}>
                              {row[c] != null ? String(row[c]) : ''}
                            </td>
                          ))}
                        </tr>
                      ))}
                    </tbody>
                    {statsData.rows.length > 1 && (
                      <tfoot>
                        <tr style={{ fontWeight: 600, borderTop: '2px solid var(--color-border)' }}>
                          {statsData.columns.map((c, ci) => {
                            const vals = statsData.rows.map(r => Number(r[c])).filter(v => !isNaN(v))
                            if (vals.length === 0) return <td key={c} style={{ padding: '4px 10px' }}>{ci === 0 ? 'Median' : ''}</td>
                            const sorted = [...vals].sort((a, b) => a - b)
                            const mid = Math.floor(sorted.length / 2)
                            const med = sorted.length % 2 !== 0 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2
                            return <td key={c} style={{ padding: '4px 10px', whiteSpace: 'nowrap' }}>{Number.isInteger(med) ? med.toLocaleString() : med.toLocaleString(undefined, { maximumFractionDigits: 1 })}</td>
                          })}
                        </tr>
                      </tfoot>
                    )}
                  </table>
                </div>
              ) : (
                <p className="loading" style={{ marginTop: 8 }}>Loading…</p>
              )
            )}
          </>
        )}

        {/* Hint */}
        <p style={{ fontSize: '.78rem', color: 'var(--color-muted-fg)', fontStyle: 'italic', margin: '8px 0 12px' }}>{hint}</p>

        {/* Config strips */}
        {configMap && configSections.map(sec => (
          <StageConfig
            key={sec.prefix}
            configMap={configMap}
            prefixes={[sec.prefix]}
            study={study}
            run={run}
            group={group}
            onConfigChanged={onConfigChanged}
          />
        ))}

        {/* Actions bar */}
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginTop: 14 }}>
          <div>
            {tabStale && <StaleKeysBadge staleKeys={allStaleKeys} configMap={configMap} />}
          </div>
          <button className="btn btn-primary" onClick={() => onRunStage(tabDef.runTarget)} disabled={isRunning}>
            {isRunning ? 'Running…' : tabStale ? `Re-run ${tabDef.label}` : `Run ${tabDef.label}`}
          </button>
        </div>
      </div>

      {/* Logs (shared across all sub-tabs) */}
      {(dada2Data?.logs.length ?? 0) > 0 && (
        <div className="card">
          <div className="card-title">Stage Logs</div>
          {dada2Data!.logs.map(log => (
            <div key={log.name} style={{ marginBottom: 4 }}>
              <button
                className="btn"
                style={{ fontSize: '.78rem', padding: '3px 10px', width: '100%', textAlign: 'left' }}
                onClick={() => fetchLog(log.url, log.name)}
              >
                {expandedLog === log.name ? 'Hide' : 'Show'} {log.name}
              </button>
              {expandedLog === log.name && logContent[log.name] && (
                <pre style={{
                  margin: '4px 0 8px', padding: 10, background: 'var(--color-surface)',
                  border: '1px solid var(--color-border)', borderRadius: 4,
                  fontSize: '.75rem', maxHeight: 300, overflow: 'auto', whiteSpace: 'pre-wrap'
                }}>
                  {logContent[log.name]}
                </pre>
              )}
            </div>
          ))}
        </div>
      )}
    </>
  )
}
