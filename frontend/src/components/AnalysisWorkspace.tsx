// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useMemo, useState, type ReactNode } from 'react'
import type { AnnotationSource, ComparisonRunSpec } from '../api/types'
import { AnalysisSourceContext } from './analysisSource'
import { useSharedResultsTables } from './annotationShared'
import { DiversityPanel } from './DiversityPanel'
import { TaxaCompositionChart } from './TaxaCompositionChart'
import { VennPanel } from './VennPanel'
import { DifferentialPanel } from './DifferentialPanel'
import { PublicationTablesPanel } from './PublicationTablesPanel'
import { FigureBuilder } from '../figure/FigureBuilder'
import { useTabParam } from '../hooks/useTabParam'
import styles from './AnalysisWorkspace.module.css'

const TABS = ['diversity', 'composition', 'overlap', 'differential', 'tables', 'figures'] as const
type AnalysisTab = (typeof TABS)[number]
const TAB_LABELS: Record<AnalysisTab, string> = {
  diversity:   'Diversity',
  composition: 'Composition',
  overlap:     'Taxon Overlap',
  differential: 'Differential Abundance',
  tables:      'Publication Tables',
  figures:     'Figures',
}

const specKey = (r: ComparisonRunSpec) => `${r.group ?? ''}|${r.run}|${r.prefix ?? ''}`

/**
 * One home for cross-run analysis: a scope bar (which runs/sub-groups, which
 * results table, whether sub-groups are pooled) shared by every tab, then tabs
 * for diversity, composition, taxon overlap, differential abundance and
 * publication tables.
 *
 * Tabs mount on first visit and stay mounted, so computed charts survive
 * switching between them.
 */
export function AnalysisWorkspace({ study, runs: baseRuns, source, perRun }: {
  study: string
  runs: ComparisonRunSpec[]
  /** Classifier whose taxonomy and categories every panel reads. */
  source: AnnotationSource
  /** Single-run views: per-sample diversity shown above the sub-group comparison,
   *  and a per-run composition chart in place of the cross-run one. */
  perRun?: { diversity?: ReactNode; composition?: ReactNode }
}) {
  const runs = useMemo(() => baseRuns.map(r => ({ ...r, source })), [baseRuns, source])
  const [tab, setTab] = useTabParam<AnalysisTab>('analysis', TABS, 'diversity')
  const [visited, setVisited] = useState<Set<AnalysisTab>>(() => new Set([tab]))
  useEffect(() => { setVisited(v => v.has(tab) ? v : new Set(v).add(tab)) }, [tab])

  //## Scope: which specs take part (all by default, including ones added later)
  const [excluded, setExcluded] = useState<Set<string>>(new Set())
  const selectedRuns = useMemo(
    () => runs.filter(r => !excluded.has(specKey(r))),
    // Keyed on content so a new-but-equal `runs` array does not refire effects.
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [runs.map(specKey).join(','), excluded],
  )
  const toggleSpec = (r: ComparisonRunSpec) => setExcluded(prev => {
    const next = new Set(prev)
    const k = specKey(r)
    next.has(k) ? next.delete(k) : next.add(k)
    return next
  })

  const [aggregate, setAggregate] = useState(false)
  const hasSubgroups = runs.some(r => r.prefix)

  // When aggregating, collapse per-sub-group specs to one per (group, run) with
  // no prefix; the backend then pools each whole run.
  const effectiveRuns = useMemo(() => {
    if (!aggregate) return selectedRuns
    const seen = new Map<string, ComparisonRunSpec>()
    for (const r of selectedRuns) {
      const key = `${r.group ?? ''}|${r.run}`
      seen.has(key) || seen.set(key, { ...r, prefix: null })
    }
    return [...seen.values()]
  }, [aggregate, selectedRuns])

  //## Dataset: results tables shared by every selected run
  const options = useSharedResultsTables(study, selectedRuns)
  const [analysisKey, setAnalysisKey] = useState<string | null>(null)
  useEffect(() => {
    setAnalysisKey(current => current && options.some(o => o.key === current) ? current : options[0]?.key ?? null)
  }, [options])
  const selected = options.find(o => o.key === analysisKey) ?? null

  // Group specs by run so sub-groups read as chips under their run.
  const byRun = useMemo(() => {
    const m = new Map<string, ComparisonRunSpec[]>()
    for (const r of runs) {
      const k = r.group ? `${r.group}/${r.run}` : r.run
      m.set(k, [...(m.get(k) ?? []), r])
    }
    return [...m.entries()]
  }, [runs])

  const usesDataset = tab === 'diversity' || tab === 'overlap' || tab === 'differential'
  const usesScope = tab !== 'figures'
  const usesAggregate = tab === 'diversity' || tab === 'composition' || tab === 'differential'

  return (
    <AnalysisSourceContext.Provider value={source}>
      <div>
        <div className="card" hidden={!usesScope}>
          <div className={styles.scope}>
            <div className={styles.row}>
              <span className={styles.rowLabel}>Compare</span>
              {byRun.map(([label, specs]) => (
                <div key={label} className={styles.chips}>
                  {specs.map(r => {
                    const on = !excluded.has(specKey(r))
                    return (
                      <button key={specKey(r)} type="button"
                        className={`${styles.chip} ${on ? styles.on : ''}`}
                        aria-pressed={on}
                        onClick={() => toggleSpec(r)}>
                        {label}{r.prefix && <span className={styles.sub}>· {r.prefix.replace(/_/g, ' ')}</span>}
                      </button>
                    )
                  })}
                </div>
              ))}
              {excluded.size > 0 && (
                <button type="button" className={styles.link} onClick={() => setExcluded(new Set())}>
                  Select all
                </button>
              )}
            </div>
            <div className={styles.row}>
              <span className={styles.rowLabel}>Dataset</span>
              {options.length > 0 ? (
                <select className={styles.select} value={analysisKey ?? ''}
                  disabled={!usesDataset}
                  title={usesDataset ? undefined : 'This tab chooses its own tables'}
                  onChange={e => setAnalysisKey(e.target.value)}>
                  {options.map(o => <option key={o.key} value={o.key}>{o.label}</option>)}
                </select>
              ) : (
                <span style={{ color: 'var(--color-muted-fg)' }}>no results table shared by the selected runs</span>
              )}
              {hasSubgroups && (
                <label className={styles.inline}
                  title={usesAggregate ? undefined : 'Not used by this tab'}
                  style={usesAggregate ? undefined : { opacity: .5 }}>
                  <input type="checkbox" checked={aggregate} disabled={!usesAggregate}
                    onChange={e => setAggregate(e.target.checked)} />
                  Aggregate sub-groups
                </label>
              )}
            </div>
          </div>
        </div>

        <div className="tabs" role="tablist">
          {TABS.map(t => (
            <button key={t} role="tab" aria-selected={tab === t} className={`tab ${tab === t ? 'active' : ''}`} onClick={() => setTab(t)}>
              {TAB_LABELS[t]}
            </button>
          ))}
        </div>

        {visited.has('diversity') && (
          <div hidden={tab !== 'diversity'}>
            {perRun?.diversity && (
              <section className={styles.section}>
                <h3 className={styles.sectionTitle}>Per sample</h3>
                {perRun.diversity}
              </section>
            )}
            {(!perRun?.diversity || selectedRuns.length >= 2) && (
              <section className={styles.section}>
                {perRun?.diversity && <h3 className={styles.sectionTitle}>Between sub-groups</h3>}
                <DiversityPanel study={study} runs={effectiveRuns} option={selected} aggregate={aggregate} />
              </section>
            )}
          </div>
        )}
        {visited.has('composition') && (
          <div hidden={tab !== 'composition'}>
            {perRun?.composition ?? (effectiveRuns.length > 0
              ? <TaxaCompositionChart study={study} runs={effectiveRuns} defaultTag="rank" />
              : <p className="empty-state">Select at least one run or sub-group.</p>)}
          </div>
        )}
        {visited.has('overlap') && (
          <div hidden={tab !== 'overlap'}>
            {selectedRuns.length >= 2
              ? <VennPanel study={study} runs={selectedRuns} option={selected} />
              : <p className="empty-state">Select at least two runs or sub-groups to compare taxon overlap.</p>}
          </div>
        )}
        {visited.has('differential') && (
          <div hidden={tab !== 'differential'}>
            {effectiveRuns.length === 2
              ? <DifferentialPanel study={study} runs={effectiveRuns} option={selected} aggregate={aggregate} />
              : <p className="empty-state">Select exactly two runs or sub-groups to test differential abundance; the first is the reference.</p>}
          </div>
        )}
        {visited.has('tables') && (
          <div hidden={tab !== 'tables'}>
            {selectedRuns.length > 0
              ? <PublicationTablesPanel study={study} runs={selectedRuns} />
              : <p className="empty-state">Select at least one run or sub-group.</p>}
          </div>
        )}
        {visited.has('figures') && (
          <div hidden={tab !== 'figures'}>
            <FigureBuilder study={study} runs={runs} source={source} />
          </div>
        )}
      </div>
    </AnalysisSourceContext.Provider>
  )
}
