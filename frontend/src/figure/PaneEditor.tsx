// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

import { useEffect, useMemo, useState } from 'react'
import { api } from '../api/client'
import { FACET_DIMENSIONS, type AnnotationSource, type CategorySet, type ComparisonRunSpec, type FacetDimension } from '../api/types'
import { errorMessage } from '../api/errorMessage'
import { useSharedResultsTables } from '../components/annotationShared'
import { useToast } from '../components/Toast'
import type { AlphaMetricKey, ChartKind, ChartSpec, FigurePane, PaneStyle, PaneText, TextLine } from './types'
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

const newText = (): PaneText => ({ lines: [{ text: '' }], align: 'center', valign: 'middle' })

export function PaneEditor({ study, figureId, runs, source, pane, onChange }: {
  study: string
  figureId: string
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
        <span>Shows</span>
        <select value={chart?.kind ?? (pane.image ? 'image' : pane.text ? 'text' : '')} onChange={e => {
          const v = e.target.value
          // A picture pane's text is its caption, so it is kept while the pane stays a picture.
          if (v === 'image') onChange({ ...pane, chart: null, image: pane.image ?? { file: '', fit: 'contain' }, text: pane.image ? pane.text : null })
          else if (v === 'text') onChange({ ...pane, chart: null, image: null, text: pane.text ?? newText() })
          else onChange({ ...pane, image: null, text: null,
                          chart: v ? defaultChart(v as ChartKind, chart, runs, source) : null })
        }}>
          <option value="">Empty</option>
          {(Object.keys(KIND_LABELS) as ChartKind[]).map(k => <option key={k} value={k}>{KIND_LABELS[k]}</option>)}
          <option value="image">Picture</option>
          <option value="text">Text</option>
        </select>
      </label>
      {!chart && pane.image && (
        <ImageFields study={study} figureId={figureId} pane={pane} onChange={onChange} />
      )}
      {!chart && pane.image && (pane.text
        ? <>
            <div className={styles.inline}>
              <span style={{ fontWeight: 600 }}>Caption</span>
              <button className="btn btn-sm" onClick={() => onChange({ ...pane, text: null })}>Remove caption</button>
            </div>
            <TextFields text={pane.text} onChange={text => onChange({ ...pane, text })} />
          </>
        : <button className="btn btn-sm" onClick={() => onChange({ ...pane, text: { ...newText(), valign: 'top' } })}>Add caption</button>)}
      {!chart && !pane.image && pane.text && <TextFields text={pane.text} onChange={text => onChange({ ...pane, text })} />}
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
                  <option value="">One pooled bar per run</option>
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
          <label className={styles.field}>
            <span>Y-axis title</span>
            <input type="text" value={pane.yTitle ?? ''} placeholder="Chart's own" disabled={pane.yTitle === ''}
              onChange={e => onChange({ ...pane, yTitle: e.target.value === '' ? null : e.target.value })} />
          </label>
          <label className={styles.inline}>
            <input type="checkbox" checked={pane.yTitle === ''}
              onChange={e => onChange({ ...pane, yTitle: e.target.checked ? '' : null })} />
            No y-axis title
          </label>
        </>
      )}
    </div>
  )
}

/** Upload a picture and choose how it fills the pane. */
function ImageFields({ study, figureId, pane, onChange }: {
  study: string; figureId: string; pane: FigurePane; onChange: (pane: FigurePane) => void
}) {
  const toast = useToast()
  const [busy, setBusy] = useState(false)
  const upload = async (f: File | undefined) => {
    if (!f) return
    setBusy(true)
    try {
      const { file } = await api.figures.uploadImage(study, figureId, f)
      // Replacing the picture keeps its caption.
      onChange({ ...pane, chart: null, image: { file, fit: pane.image?.fit ?? 'contain' } })
    } catch (e) {
      toast.error(`Picture not uploaded: ${errorMessage(e)}`)
    } finally {
      setBusy(false)
    }
  }
  return (
    <>
      <label className={styles.field}>
        <span>Picture</span>
        <input type="file" accept="image/png,image/jpeg" disabled={busy}
          onChange={e => { upload(e.target.files?.[0]); e.target.value = '' }} />
      </label>
      {busy && <span style={{ color: 'var(--color-muted-fg)' }}>Uploading…</span>}
      {pane.image?.file && (
        <label className={styles.field}>
          <span>Fit</span>
          <select value={pane.image.fit} onChange={e => onChange({ ...pane, image: { ...pane.image!, fit: e.target.value as 'contain' | 'cover' } })}>
            <option value="contain">Whole picture</option>
            <option value="cover">Fill the pane, cropping</option>
          </select>
        </label>
      )}
      <span style={{ color: 'var(--color-muted-fg)' }}>PNG or JPEG. A PNG keeps its transparency.</span>
    </>
  )
}

/** Lines of text, each with its own size and colour. */
function TextFields({ text, onChange }: { text: PaneText; onChange: (text: PaneText) => void }) {
  const setLine = (i: number, patch: Partial<TextLine>) =>
    onChange({ ...text, lines: text.lines.map((l, k) => k === i ? { ...l, ...patch } : l) })
  return (
    <>
      <div className={styles.field}>
        <span>Align</span>
        <span className={styles.inline}>
          <select value={text.align} aria-label="Horizontal alignment"
            onChange={e => onChange({ ...text, align: e.target.value as PaneText['align'] })}>
            <option value="left">Left</option>
            <option value="center">Centre</option>
            <option value="right">Right</option>
          </select>
          <select value={text.valign} aria-label="Vertical alignment"
            onChange={e => onChange({ ...text, valign: e.target.value as PaneText['valign'] })}>
            <option value="top">Top</option>
            <option value="middle">Middle</option>
            <option value="bottom">Bottom</option>
          </select>
          <select value={text.rotate ?? 0} aria-label="Rotation"
            onChange={e => onChange({ ...text, rotate: Number(e.target.value) === 0 ? undefined : Number(e.target.value) as 90 | -90 })}>
            <option value="0">Upright</option>
            <option value="-90">Turned, reading up</option>
            <option value="90">Turned, reading down</option>
          </select>
        </span>
      </div>
      <div className={styles.field}>
        <span>Move (mm)</span>
        <span className={styles.inline}>
          {(['right', 'down'] as const).map((dir, k) => (
            <input key={dir} type="number" step={0.1} style={{ width: 60 }} aria-label={`Move ${dir}, in mm`} placeholder={dir}
              value={text.offsetMm?.[k] ?? ''}
              onChange={e => {
                const o: [number, number] = [...(text.offsetMm ?? [0, 0])] as [number, number]
                o[k] = Number(e.target.value) || 0
                onChange({ ...text, offsetMm: o[0] || o[1] ? o : undefined })
              }} />
          ))}
        </span>
      </div>
      {text.lines.map((l, i) => (
        <div key={i} className={styles.inline} style={{ flexWrap: 'nowrap' }}>
          <input type="text" value={l.text} aria-label={`Line ${i + 1}`} style={{ flex: 1, minWidth: 0, font: 'inherit' }}
            onChange={e => setLine(i, { text: e.target.value })} />
          <input type="number" min={4} max={72} step={0.5} style={{ width: 52 }} aria-label={`Size of line ${i + 1} in pt`}
            value={l.sizePt ?? ''} placeholder="pt"
            onChange={e => setLine(i, { sizePt: e.target.value === '' ? undefined : Math.min(72, Math.max(4, Number(e.target.value) || 8)) })} />
          <input type="color" value={l.colour ?? '#222222'} aria-label={`Colour of line ${i + 1}`}
            onChange={e => setLine(i, { colour: e.target.value })} />
          <button className="btn btn-sm" aria-label={`Remove line ${i + 1}`} disabled={text.lines.length === 1}
            onClick={() => onChange({ ...text, lines: text.lines.filter((_, k) => k !== i) })}>✕</button>
        </div>
      ))}
      <div className={styles.inline}>
        <button className="btn btn-sm" onClick={() => onChange({ ...text, lines: [...text.lines, { text: '' }] })}>Add line</button>
        <span style={{ color: 'var(--color-muted-fg)' }}>*italic*, **bold**. Blank size uses the text size.</span>
      </div>
    </>
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
      <div className={styles.inline}>
        <label className={styles.inline}>
          <input type="checkbox" checked={style.hideXTicks ?? false} onChange={e => set('hideXTicks', e.target.checked || undefined)} />
          No x tick labels
        </label>
        <label className={styles.inline}>
          <input type="checkbox" checked={style.hideYTicks ?? false} onChange={e => set('hideYTicks', e.target.checked || undefined)} />
          No y tick labels
        </label>
      </div>
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
      {chart?.kind === 'alpha' && (
        <label className={styles.field}>
          <span>Box fill opacity</span>
          <input type="number" min={0} max={1} step={0.05} value={style.fillOpacity ?? ''} placeholder="Chart's own"
            onChange={e => set('fillOpacity', numberOrUndefined(e.target.value, 0, 1))} />
        </label>
      )}
      {hasPoints && (
        <label className={styles.field}>
          <span>Point opacity</span>
          <input type="number" min={0} max={1} step={0.05} value={style.pointOpacity ?? ''} placeholder="Chart's own"
            onChange={e => set('pointOpacity', numberOrUndefined(e.target.value, 0, 1))} />
        </label>
      )}
      {(chart?.kind === 'composition' || chart?.kind === 'alpha') && (
        <label className={styles.inline}>
          <input type="checkbox" checked={style.horizontal ?? false} onChange={e => set('horizontal', e.target.checked || undefined)} />
          Horizontal
        </label>
      )}
      {chart?.kind === 'composition' && (
        <label className={styles.inline}>
          <input type="checkbox" checked={style.labelSegments ?? false} onChange={e => set('labelSegments', e.target.checked || undefined)} />
          Name segments inside the bars
        </label>
      )}
      <label className={styles.field}>
        <span>Comparison</span>
        <input type="text" value={style.compareLabel ?? ''} placeholder="e.g. p = 0.015"
          onChange={e => set('compareLabel', e.target.value || undefined)} />
      </label>
      <label className={styles.inline}>
        <input type="checkbox" checked={style.hideFromLegend ?? false} onChange={e => set('hideFromLegend', e.target.checked || undefined)} />
        Leave out of the group legend
      </label>
      <label className={styles.inline}>
        <input type="checkbox" checked={style.categoryBands ?? false} onChange={e => set('categoryBands', e.target.checked || undefined)} />
        Band each category in its colour
      </label>
      <label className={styles.field}>
        <span>Gridlines</span>
        <select value={style.grid ?? ''} onChange={e => set('grid', (e.target.value || undefined) as PaneStyle['grid'])}>
          <option value="">Chart's own</option>
          <option value="strong">Darker</option>
          <option value="none">None</option>
        </select>
      </label>
      {style.valueRange && (
        <label className={styles.field}>
          <span>Tick step</span>
          <input type="number" step="any" min={0} value={style.valueStep ?? ''} placeholder="ends and middle"
            onChange={e => set('valueStep', numberOrUndefined(e.target.value, 0, 1e9) || undefined)} />
        </label>
      )}
      <div className={styles.field}>
        <span>Value axis</span>
        <span className={styles.inline}>
          <input type="number" step="any" value={style.valueRange?.[0] ?? ''} placeholder="from" aria-label="Value axis start" style={{ width: 60 }}
            onChange={e => { const v = numberOrUndefined(e.target.value, -1e9, 1e9); set('valueRange', v == null ? undefined : [v, style.valueRange?.[1] ?? v + 1]) }} />
          <input type="number" step="any" value={style.valueRange?.[1] ?? ''} placeholder="to" aria-label="Value axis end" style={{ width: 60 }}
            onChange={e => { const v = numberOrUndefined(e.target.value, -1e9, 1e9); set('valueRange', v == null ? undefined : [style.valueRange?.[0] ?? 0, v]) }} />
        </span>
      </div>
      <label className={styles.inline}>
        <input type="checkbox" checked={style.hideNotes ?? false} onChange={e => set('hideNotes', e.target.checked || undefined)} />
        No chart notes (test results, warnings)
      </label>
      {isGrid && (
        <label className={styles.inline}>
          <input type="checkbox" checked={style.hideRowHeaders ?? false} onChange={e => set('hideRowHeaders', e.target.checked || undefined)} />
          No row headers
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
