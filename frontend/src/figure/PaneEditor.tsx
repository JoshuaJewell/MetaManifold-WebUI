// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useMemo, useState } from 'react'
import { api } from '../api/client'
import { FACET_DIMENSIONS, type AnnotationSource, type CategorySet, type ComparisonRunSpec, type FacetDimension } from '../api/types'
import { useSharedResultsTables } from '../components/annotationShared'
import type { AlphaMetricKey, ChartKind, ChartSpec, FigurePane, PaneStyle } from './types'
import styles from './Figure.module.css'

const specKey = (r: ComparisonRunSpec) => `${r.group ?? ''}|${r.run}|${r.prefix ?? ''}`
const specLabel = (r: ComparisonRunSpec) =>
  `${r.group ? `${r.group}/` : ''}${r.run}${r.prefix ? ` · ${r.prefix.replace(/_/g, ' ')}` : ''}`

const KIND_LABELS: Record<ChartKind, string> = {
  alpha: 'Alpha diversity',
  nmds: 'NMDS',
  composition: 'Composition bars',
}

function defaultChart(kind: ChartKind, prev: ChartSpec | null, runs: ComparisonRunSpec[], source: AnnotationSource): ChartSpec {
  const base = {
    source: prev?.source ?? source,
    runs: prev?.runs ?? runs.map(({ run, group, prefix }) => ({ run, group, prefix })),
    table: prev?.table ?? 'merged',
    aggregate: prev?.aggregate ?? false,
  }
  if (kind === 'alpha') return { ...base, kind, metric: 'richness' }
  if (kind === 'nmds') return { ...base, kind }
  return { ...base, kind, tag: 'rank', value: 'Genus', topN: 15, relative: true, mode: 'stacked', keepEmpty: false, facet: null }
}

export function PaneEditor({ study, runs, source, pane, onChange }: {
  study: string
  runs: ComparisonRunSpec[]
  source: AnnotationSource
  pane: FigurePane
  onChange: (pane: FigurePane) => void
}) {
  const chart = pane.chart
  const set = (patch: Partial<ChartSpec>) => chart && onChange({ ...pane, chart: { ...chart, ...patch } as ChartSpec })

  const selected = useMemo(() => {
    const keys = new Set((chart?.runs ?? []).map(specKey))
    return runs.filter(r => keys.has(specKey(r)))
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [runs, chart?.runs.map(specKey).join(',')])
  const tables = useSharedResultsTables(study, selected)

  const [ranks, setRanks] = useState<string[]>([])
  const [catSets, setCatSets] = useState<CategorySet[]>([])
  const first = chart?.runs[0]
  const isComposition = chart?.kind === 'composition'
  const tag = chart?.kind === 'composition' ? chart.tag : null
  useEffect(() => {
    if (!isComposition || !first || !chart) return
    if (tag === 'rank') {
      api.analysis.ranks(study, first.run, { table: chart.table, group: first.group ?? undefined, source: chart.source })
        .then(setRanks).catch(() => setRanks([]))
    } else {
      api.composition.categorySets().then(setCatSets).catch(() => setCatSets([]))
    }
  }, [study, isComposition, tag, first?.run, first?.group, chart?.table, chart?.source]) // eslint-disable-line react-hooks/exhaustive-deps

  const toggleRun = (r: ComparisonRunSpec) => {
    if (!chart) return
    const on = chart.runs.some(x => specKey(x) === specKey(r))
    const next = on ? chart.runs.filter(x => specKey(x) !== specKey(r))
                    : runs.filter(x => specKey(x) === specKey(r) || chart.runs.some(c => specKey(c) === specKey(x)))
                          .map(({ run, group, prefix }) => ({ run, group, prefix }))
    set({ runs: next })
  }

  return (
    <div className={styles.section}>
      <div className={styles.sectionTitle}>Pane</div>
      <label className={styles.field}>
        <span>Chart</span>
        <select value={chart?.kind ?? ''} onChange={e => onChange({
          ...pane, chart: e.target.value ? defaultChart(e.target.value as ChartKind, chart, runs, source) : null,
        })}>
          <option value="">Empty</option>
          {(Object.keys(KIND_LABELS) as ChartKind[]).map(k => <option key={k} value={k}>{KIND_LABELS[k]}</option>)}
        </select>
      </label>
      {chart && (
        <>
          <label className={styles.field}>
            <span>Classifier</span>
            <select value={chart.source} onChange={e => set({ source: e.target.value as AnnotationSource })}>
              <option value="DADA2">DADA2</option>
              <option value="VSEARCH">VSEARCH</option>
            </select>
          </label>
          <div className={styles.field} style={{ alignItems: 'start' }}>
            <span>Compare</span>
            <div className={styles.chips}>
              {runs.map(r => {
                const on = chart.runs.some(x => specKey(x) === specKey(r))
                return (
                  <button key={specKey(r)} type="button" aria-pressed={on}
                    className={`${styles.chip} ${on ? styles.on : ''}`} onClick={() => toggleRun(r)}>
                    {specLabel(r)}
                  </button>
                )
              })}
            </div>
          </div>
          <label className={styles.field}>
            <span>Table</span>
            <select value={chart.table} onChange={e => set({ table: e.target.value })}>
              {!tables.some(t => t.table === chart.table) && <option value={chart.table}>{chart.table}</option>}
              {tables.map(t => <option key={t.key} value={t.table}>{t.label}</option>)}
            </select>
          </label>
          {runs.some(r => r.prefix) && chart.kind !== 'composition' && (
            <label className={styles.inline}>
              <input type="checkbox" checked={chart.aggregate} onChange={e => set({ aggregate: e.target.checked })} />
              Pool each run's sub-groups
            </label>
          )}

          {chart.kind === 'alpha' && (
            <label className={styles.field}>
              <span>Metric</span>
              <select value={chart.metric} onChange={e => set({ metric: e.target.value as AlphaMetricKey })}>
                <option value="richness">Richness</option>
                <option value="shannon">Shannon</option>
                <option value="simpson">Simpson</option>
              </select>
            </label>
          )}

          {chart.kind === 'composition' && (
            <>
              <label className={styles.field}>
                <span>Bars by</span>
                <select value={chart.tag} onChange={e => set({ tag: e.target.value as 'rank' | 'category',
                                                             value: e.target.value === 'rank' ? 'Genus' : 'default' })}>
                  <option value="rank">Taxonomic rank</option>
                  <option value="category">Category set</option>
                </select>
              </label>
              <label className={styles.field}>
                <span>{chart.tag === 'rank' ? 'Rank' : 'Set'}</span>
                <select value={chart.value} onChange={e => set({ value: e.target.value })}>
                  {chart.tag === 'rank'
                    ? (ranks.includes(chart.value) ? ranks : [chart.value, ...ranks]).map(r => <option key={r} value={r}>{r}</option>)
                    : (catSets.some(c => c.name === chart.value) ? catSets : [{ name: chart.value, label: chart.value } as CategorySet, ...catSets])
                        .map(c => <option key={c.name} value={c.name}>{c.label}</option>)}
                </select>
              </label>
              {chart.tag === 'rank' && (
                <label className={styles.field}>
                  <span>Top taxa</span>
                  <input type="number" min={1} max={50} value={chart.topN}
                    onChange={e => set({ topN: Math.max(1, Number(e.target.value) || 1) })} />
                </label>
              )}
              <label className={styles.field}>
                <span>Grid</span>
                <select value={chart.facet ? `${chart.facet.rows}|${chart.facet.cols}` : ''}
                  onChange={e => {
                    const [rows, cols] = e.target.value.split('|')
                    set({ facet: e.target.value ? { rows: rows as FacetDimension, cols: cols as FacetDimension } : null })
                  }}>
                  <option value="">One chart</option>
                  {FACET_DIMENSIONS.flatMap(r => FACET_DIMENSIONS.filter(c => c !== r).map(c => (
                    <option key={`${r}|${c}`} value={`${r}|${c}`}>Rows by {r}, columns by {c}</option>
                  )))}
                </select>
              </label>
              <div className={styles.inline}>
                <label className={styles.inline}>
                  <input type="checkbox" checked={chart.relative} onChange={e => set({ relative: e.target.checked })} />
                  Relative
                </label>
                <label className={styles.inline}>
                  <input type="checkbox" checked={chart.mode === 'grouped'}
                    onChange={e => set({ mode: e.target.checked ? 'grouped' : 'stacked' })} />
                  Grouped bars
                </label>
                <label className={styles.inline}>
                  <input type="checkbox" checked={chart.keepEmpty} onChange={e => set({ keepEmpty: e.target.checked })} />
                  Show empty as gaps
                </label>
              </div>
            </>
          )}

          <label className={styles.field}>
            <span>Title</span>
            <input type="text" value={pane.title ?? ''} placeholder="Chart's own title" disabled={pane.title === ''}
              onChange={e => onChange({ ...pane, title: e.target.value === '' ? null : e.target.value })} />
          </label>
          <label className={styles.inline}>
            <input type="checkbox" checked={pane.title === ''}
              onChange={e => onChange({ ...pane, title: e.target.checked ? '' : null })} />
            No title
          </label>
          <PaneStyleFields pane={pane} onChange={onChange} />
          <label className={styles.field}>
            <span>X-axis title</span>
            <input type="text" value={pane.xTitle ?? ''} placeholder="Chart's own"
              onChange={e => onChange({ ...pane, xTitle: e.target.value === '' ? null : e.target.value })} />
          </label>
        </>
      )}
    </div>
  )
}

/** Optional adjustments to the chart's look; blank keeps the chart's own. */
function PaneStyleFields({ pane, onChange }: { pane: FigurePane; onChange: (pane: FigurePane) => void }) {
  const style = pane.style ?? {}
  const set = <K extends keyof PaneStyle>(key: K, value: PaneStyle[K] | undefined) => {
    const next = { ...style }
    if (value === undefined) delete next[key]
    else next[key] = value
    onChange({ ...pane, style: next })
  }
  const numberOrUndefined = (v: string, min: number, max: number) =>
    v === '' || !Number.isFinite(Number(v)) ? undefined : Math.min(max, Math.max(min, Number(v)))
  const chart = pane.chart
  const hasLines = chart?.kind === 'alpha'
  const hasPoints = chart?.kind === 'alpha' || chart?.kind === 'nmds'
  const isGrid = chart?.kind === 'composition' && chart.facet != null
  return (
    <>
      <label className={styles.field}>
        <span>X tick angle</span>
        <select value={style.tickAngle ?? ''} onChange={e => set('tickAngle', e.target.value === '' ? undefined : Number(e.target.value))}>
          <option value="">Chart's own</option>
          <option value="0">Horizontal</option>
          <option value="-45">45°</option>
          <option value="-90">Vertical</option>
        </select>
      </label>
      <label className={styles.inline}>
        <input type="checkbox" checked={style.allTicks ?? false} onChange={e => set('allTicks', e.target.checked || undefined)} />
        Label every x category
      </label>
      {hasLines && (
        <label className={styles.field}>
          <span>Pair lines opacity</span>
          <input type="number" min={0} max={1} step={0.05} value={style.lineOpacity ?? ''} placeholder="Chart's own"
            onChange={e => set('lineOpacity', numberOrUndefined(e.target.value, 0, 1))} />
        </label>
      )}
      {hasPoints && (
        <label className={styles.field}>
          <span>Point outline (px)</span>
          <input type="number" min={0} max={4} step={0.1} value={style.pointBorder ?? ''} placeholder="Chart's own"
            onChange={e => set('pointBorder', numberOrUndefined(e.target.value, 0, 4))} />
        </label>
      )}
      {hasPoints && (
        <label className={styles.field}>
          <span>Point opacity</span>
          <input type="number" min={0} max={1} step={0.05} value={style.pointOpacity ?? ''} placeholder="Chart's own"
            onChange={e => set('pointOpacity', numberOrUndefined(e.target.value, 0, 1))} />
        </label>
      )}
      {isGrid && (
        <label className={styles.field}>
          <span>Column spacing</span>
          <input type="number" min={0} max={0.4} step={0.01} value={style.gridGap ?? ''} placeholder="Chart's own"
            onChange={e => set('gridGap', numberOrUndefined(e.target.value, 0, 0.4))} />
        </label>
      )}
      {isGrid && (
        <label className={styles.field}>
          <span>Row spacing</span>
          <input type="number" min={0} max={0.4} step={0.01} value={style.gridRowGap ?? ''} placeholder="Chart's own"
            onChange={e => set('gridRowGap', numberOrUndefined(e.target.value, 0, 0.4))} />
        </label>
      )}
    </>
  )
}
