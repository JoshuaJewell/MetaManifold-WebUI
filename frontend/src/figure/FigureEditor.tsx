// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useLayoutEffect, useRef, useState } from 'react'
import Plotly from 'plotly.js-dist-min'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import type { AnnotationSource, ComparisonRunSpec } from '../api/types'
import { AddToReport } from '../components/AddToReport'
import { useToast } from '../components/Toast'
import { saveBlob } from '../utils/download'
import { clearChartCache, fetchChart, paneFigure, parseColour, type PlotSpec } from './charts'
import { figurePdf, figurePng, figureSvg, figureTiff, imageDataUrl } from './export'
import { layoutLegend, mergeEntries, paneLegendEntries, textWidth, type LegendEntry } from './legend'
import { PaneEditor } from './PaneEditor'
import {
  PAGE_SIZES, PX_PER_MM, PX_PER_PT, captionSplit, figureGeometry, layoutText, textFrame, legendPlacement, newGroup, resizeGroup, tickPt, titlePt,
  type FigureDoc, type FigureGroup, type FigurePane, type LegendPlacement, type PageSize,
} from './types'
import styles from './Figure.module.css'

/** Index of the pane that carries a group's shared legend: the top-right pane with a chart. */
function legendPane(g: FigureGroup): number {
  for (let r = 0; r < g.rows; r++)
    for (let c = g.cols - 1; c >= 0; c--)
      if (g.panes[r * g.cols + c]?.chart) return r * g.cols + c
  return -1
}

const showsLegend = (g: FigureGroup, i: number) => {
  const place = legendPlacement(g)
  return place === 'each' || (place === 'right' && legendPane(g) === i)
}

const fontFamily = (doc: FigureDoc) =>
  doc.font.family === 'Arial' ? 'Arial, Helvetica, sans-serif' : '"Times New Roman", Times, serif'

/** A group's horizontal legend, positioned inside its strip. */
function GroupLegend({ entries, doc, box, italic, bold }: {
  entries: LegendEntry[]; doc: FigureDoc; box: { x: number; y: number; w: number; h: number }; italic?: boolean; bold?: boolean
}) {
  const family = fontFamily(doc)
  const w = box.w * PX_PER_MM, h = box.h * PX_PER_MM
  const lay = layoutLegend(entries, w, doc.font.sizePt * PX_PER_PT, family, italic, bold)
  return (
    <div style={{ position: 'absolute', left: box.x * PX_PER_MM, top: box.y * PX_PER_MM, width: w, height: h,
                  fontFamily: family, fontSize: lay.fontPx, fontStyle: italic ? 'italic' : undefined,
                  fontWeight: bold ? 'bold' : undefined, color: '#444', pointerEvents: 'none' }}>
      {lay.items.map(it => (
        <span key={it.entry.name}>
          <span style={{ position: 'absolute', left: it.x, top: (h - lay.swatch) / 2, width: lay.swatch, height: lay.swatch,
                         background: it.entry.fill, boxSizing: 'border-box',
                         border: it.entry.stroke ? `1px solid ${it.entry.stroke}` : undefined,
                         borderRadius: it.entry.shape === 'circle' ? '50%' : 0 }} />
          <span style={{ position: 'absolute', left: it.textX, top: 0, lineHeight: `${h}px`, whiteSpace: 'nowrap',
                         fontStyle: italic || it.entry.italic ? 'italic' : undefined }}>{it.entry.name}</span>
        </span>
      ))}
    </div>
  )
}

/** Every pane's chart as it will be drawn, for export. */
async function buildPanes(study: string, doc: FigureDoc): Promise<(PlotSpec | null)[][]> {
  const geo = figureGeometry(doc)
  return Promise.all(doc.groups.map((g, gi) => Promise.all(g.panes.map(async (pane, pi) => {
    if (!pane.chart) return null
    const box = geo.groups[gi].panes[pi]
    const raw = await fetchChart(study, pane.chart)
    return paneFigure(raw, pane, doc, {
      showLegend: showsLegend(g, pi),
      widthPx: box.w * PX_PER_MM, heightPx: box.h * PX_PER_MM,
    })
  }))))
}

/** Every image pane's picture as a data URL, keyed by file. */
async function loadImages(study: string, doc: FigureDoc): Promise<Record<string, string>> {
  const files = new Set(doc.groups.flatMap(g => g.panes.flatMap(p => !p.chart && p.image?.file ? [p.image.file] : [])))
  const entries = await Promise.all([...files].map(async f =>
    [f, await imageDataUrl(api.figures.imageUrl(study, doc.id, f))] as const))
  return Object.fromEntries(entries)
}

/** The page as one SVG, with every chart fetched and every picture inlined. */
async function pageSvg(study: string, doc: FigureDoc) {
  const [panes, images] = await Promise.all([buildPanes(study, doc), loadImages(study, doc)])
  return figureSvg(doc, panes, images)
}

/** A text pane, drawn the way the export draws it. */
function PaneTextView({ pane, doc, w, h }: { pane: FigurePane; doc: FigureDoc; w: number; h: number }) {
  if (!pane.text) return null
  const frame = textFrame(pane.text, w, h)
  return (
    <svg width={w} height={h} style={{ display: 'block', overflow: 'visible' }}>
      <g transform={frame.transform}>
      {layoutText(pane.text, doc, frame.w, frame.h, (r, s) => textWidth(r.text, s, fontFamily(doc), r.italic, r.bold)).map((l, i) => (
        <text key={i} x={l.x} y={l.y} fontFamily={fontFamily(doc)} fontSize={l.sizePx} fill={l.fill} textAnchor={l.anchor}>
          {l.runs.map((r, k) => (
            <tspan key={k} fontWeight={r.bold ? 'bold' : undefined} fontStyle={r.italic ? 'italic' : undefined}>{r.text}</tspan>
          ))}
        </text>
      ))}
      </g>
    </svg>
  )
}

function PaneContent(props: {
  study: string; figureId: string; pane: FigurePane; doc: FigureDoc
  showLegend: boolean; w: number; h: number; version: number
  onEntries: (entries: LegendEntry[]) => void
}) {
  const { study, figureId, pane, doc, w, h } = props
  if (pane.chart) return <PanePlot {...props} />
  if (pane.image && !pane.image.file) return <div className={styles.paneEmpty}>Choose a picture for this pane.</div>
  if (pane.image) {
    const { imageH, captionH } = captionSplit(pane, doc, h)
    return (
      <>
        <img src={api.figures.imageUrl(study, figureId, pane.image.file)} alt="" draggable={false}
          style={{ display: 'block', width: w, height: imageH, objectFit: pane.image.fit }} />
        {captionH > 0 && <PaneTextView pane={pane} doc={doc} w={w} h={captionH} />}
      </>
    )
  }
  if (pane.text) return <PaneTextView pane={pane} doc={doc} w={w} h={h} />
  return <div className={styles.paneEmpty}>Empty pane. Select it to choose a chart, picture or text.</div>
}

function PanePlot({ study, pane, doc, showLegend, w, h, version, onEntries }: {
  study: string; pane: FigurePane; doc: FigureDoc
  showLegend: boolean; w: number; h: number; version: number
  onEntries: (entries: LegendEntry[]) => void
}) {
  const ref = useRef<HTMLDivElement>(null)
  const [error, setError] = useState<string | null>(null)
  const [loading, setLoading] = useState(false)
  const key = JSON.stringify([pane, doc.font, doc.colours, showLegend, Math.round(w), Math.round(h), version])
  useEffect(() => {
    const el = ref.current
    if (!el || !pane.chart) return
    let cancelled = false
    setLoading(true)
    setError(null)
    fetchChart(study, pane.chart)
      .then(raw => {
        if (cancelled) return
        const fig = paneFigure(raw, pane, doc, { showLegend, widthPx: w, heightPx: h })
        onEntries(paneLegendEntries(fig, pane))
        return Plotly.react(el, fig.data as Plotly.Data[], fig.layout, { staticPlot: true, displayModeBar: false })
      })
      .catch(e => { if (!cancelled) setError(errorMessage(e)) })
      .finally(() => { if (!cancelled) setLoading(false) })
    return () => { cancelled = true }
  }, [key]) // eslint-disable-line react-hooks/exhaustive-deps
  useEffect(() => () => { if (ref.current) Plotly.purge(ref.current) }, [])

  if (!pane.chart) return <div className={styles.paneEmpty}>Empty pane. Select it to choose a chart.</div>
  return (
    <>
      <div ref={ref} style={{ width: w, height: h }} />
      {loading && <div className={styles.paneEmpty} style={{ border: 'none' }}>Loading…</div>}
      {error && <div className={`${styles.paneEmpty} ${styles.paneError}`}>{error}</div>}
    </>
  )
}

type Selection = { group: number; pane: number | null }

export function FigureEditor({ study, runs, source, doc, onChange }: {
  study: string
  runs: ComparisonRunSpec[]
  source: AnnotationSource
  doc: FigureDoc
  onChange: (doc: FigureDoc) => void
}) {
  const toast = useToast()
  const [sel, setSel] = useState<Selection>({ group: 0, pane: null })
  const [version, setVersion] = useState(0)
  const dpi = doc.rasterDpi ?? 300
  const [busy, setBusy] = useState<string | null>(null)
  const [savedImage, setSavedImage] = useState<string | null>(null)
  // Legend entries each pane reported, keyed "group:pane".
  const [paneEntries, setPaneEntries] = useState<Record<string, LegendEntry[]>>({})
  const reportEntries = (key: string, entries: LegendEntry[]) => setPaneEntries(prev =>
    JSON.stringify(prev[key]) === JSON.stringify(entries) ? prev : { ...prev, [key]: entries })

  const stageRef = useRef<HTMLDivElement>(null)
  const [scale, setScale] = useState(1)
  const geo = figureGeometry(doc)
  const pageW = geo.page.w * PX_PER_MM, pageH = geo.page.h * PX_PER_MM
  useLayoutEffect(() => {
    const el = stageRef.current
    if (!el) return
    const fit = () => setScale(Math.min(1.5, (el.clientWidth - 24) / pageW))
    fit()
    const ro = new ResizeObserver(fit)
    ro.observe(el)
    return () => ro.disconnect()
  }, [pageW])

  const group = doc.groups[Math.min(sel.group, doc.groups.length - 1)]
  const setGroup = (i: number, g: FigureGroup) => onChange({ ...doc, groups: doc.groups.map((x, k) => k === i ? g : x) })
  const setPane = (gi: number, pi: number, p: FigurePane) => {
    const g = doc.groups[gi]
    setGroup(gi, { ...g, panes: g.panes.map((x, k) => k === pi ? p : x) })
  }
  const moveGroup = (i: number, step: number) => {
    const j = i + step
    if (j < 0 || j >= doc.groups.length) return
    const groups = [...doc.groups]
    ;[groups[i], groups[j]] = [groups[j], groups[i]]
    onChange({ ...doc, groups })
    setSel({ group: j, pane: null })
  }

  const exportAs = async (format: 'pdf' | 'png' | 'tif'): Promise<{ blob: Blob; ext: string }> => {
    const svg = await pageSvg(study, doc)
    const blob = format === 'pdf' ? await figurePdf(doc, svg)
               : format === 'png' ? await figurePng(doc, svg, dpi)
               : await figureTiff(doc, svg, dpi)
    return { blob, ext: `.${format}` }
  }
  const download = async (format: 'pdf' | 'png' | 'tif') => {
    setBusy(format)
    try {
      const { blob, ext } = await exportAs(format)
      saveBlob(blob, `${doc.title.replace(/[^\w.-]+/g, '_') || 'figure'}${ext}`)
    } catch (e) {
      toast.error(`Export failed: ${errorMessage(e)}`)
    } finally {
      setBusy(null)
    }
  }

  // A PNG of the page is kept beside the layout, redrawn once edits settle.
  const renderGen = useRef(0)
  useEffect(() => {
    const gen = ++renderGen.current
    const t = window.setTimeout(async () => {
      try {
        const svg = await pageSvg(study, doc)
        if (gen !== renderGen.current) return
        const png = await figurePng(doc, svg, dpi)
        if (gen !== renderGen.current) return
        const { file } = await api.figures.saveRaster(study, doc.id, png)
        setSavedImage(`${file} at ${new Date().toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}`)
      } catch (e) {
        if (gen === renderGen.current) setSavedImage(`Image not saved: ${errorMessage(e)}`)
      }
    }, 2500)
    return () => window.clearTimeout(t)
  }, [doc, version, study]) // eslint-disable-line react-hooks/exhaustive-deps

  const page = PAGE_SIZES[doc.page.size]
  /** The group's letter, or its position when it has none. */
  const groupName = (i: number) => geo.groups[i]?.letter || String(i + 1)
  const num = (v: string, min: number, max: number, fallback: number) =>
    Math.min(max, Math.max(min, Number.isFinite(Number(v)) && v !== '' ? Number(v) : fallback))

  return (
    <div>
      <div className={styles.bar}>
        <button className="btn" disabled={busy !== null} onClick={() => download('pdf')}>
          {busy === 'pdf' ? 'Exporting…' : 'PDF'}
        </button>
        <button className="btn" disabled={busy !== null} onClick={() => download('png')}>
          {busy === 'png' ? 'Exporting…' : 'PNG'}
        </button>
        <button className="btn" disabled={busy !== null} onClick={() => download('tif')}>
          {busy === 'tif' ? 'Exporting…' : 'TIFF'}
        </button>
        <label className={styles.inline}>
          at
          <select value={dpi} onChange={e => onChange({ ...doc, rasterDpi: Number(e.target.value) })} aria-label="Raster resolution">
            {[300, 600, 1200].map(d => <option key={d} value={d}>{d} dpi</option>)}
          </select>
        </label>
        <AddToReport study={study} kind="figure" className="btn" defaultTitle={doc.title}
          make={() => exportAs('pdf')} />
        <span className={styles.spacer} />
        {savedImage && <span style={{ color: 'var(--color-muted-fg)', fontSize: '.78rem' }} title="Saved in the study's figures folder">{savedImage}</span>}
        <button className="btn" title="Fetch every chart again"
          onClick={() => { clearChartCache(); setVersion(v => v + 1) }}>Refresh data</button>
      </div>

      <div className={styles.editor}>
        <div className={styles.side}>
          <div className={styles.section}>
            <div className={styles.sectionTitle}>Page</div>
            <label className={styles.field}>
              <span>Size</span>
              <select value={doc.page.size} onChange={e => {
                const size = e.target.value as PageSize
                const p = PAGE_SIZES[size]
                onChange({ ...doc, page: { ...doc.page, size,
                  heightMm: p.fixed ? doc.page.heightMm : p.h,
                  marginMm: p.fixed ? Math.max(doc.page.marginMm, 10) : 0 } })
              }}>
                {(Object.keys(PAGE_SIZES) as PageSize[]).map(k => <option key={k} value={k}>{PAGE_SIZES[k].label}</option>)}
              </select>
            </label>
            {doc.page.size === 'custom' && (
              <label className={styles.field}>
                <span>Width (mm)</span>
                <input type="number" min={30} max={600} value={doc.page.widthMm}
                  onChange={e => onChange({ ...doc, page: { ...doc.page, widthMm: num(e.target.value, 30, 600, 180) } })} />
              </label>
            )}
            {!page.fixed && (
              <label className={styles.field}>
                <span>Height (mm)</span>
                <input type="number" min={20} max={600} value={doc.page.heightMm}
                  onChange={e => onChange({ ...doc, page: { ...doc.page, heightMm: num(e.target.value, 20, 600, 200) } })} />
              </label>
            )}
            <label className={styles.field}>
              <span>Margin (mm)</span>
              <input type="number" min={0} max={40} value={doc.page.marginMm}
                onChange={e => onChange({ ...doc, page: { ...doc.page, marginMm: num(e.target.value, 0, 40, 0) } })} />
            </label>
            <label className={styles.field}>
              <span>Font</span>
              <select value={doc.font.family}
                onChange={e => onChange({ ...doc, font: { ...doc.font, family: e.target.value as FigureDoc['font']['family'] } })}>
                <option value="Arial">Arial / Helvetica</option>
                <option value="Times New Roman">Times</option>
              </select>
            </label>
            <label className={styles.field}>
              <span>Text size (pt)</span>
              <input type="number" min={4} max={24} step={0.5} value={doc.font.sizePt}
                onChange={e => onChange({ ...doc, font: { ...doc.font, sizePt: num(e.target.value, 4, 24, 8) } })} />
            </label>
            <label className={styles.field}>
              <span>Titles (pt)</span>
              <input type="number" min={4} max={36} step={0.5} value={titlePt(doc)}
                onChange={e => onChange({ ...doc, font: { ...doc.font, titlePt: num(e.target.value, 4, 36, 10) } })} />
            </label>
            <label className={styles.field}>
              <span>Bold titles</span>
              <input type="checkbox" checked={doc.font.titleBold ?? false}
                onChange={e => onChange({ ...doc, font: { ...doc.font, titleBold: e.target.checked || undefined } })} />
            </label>
            <label className={styles.field}>
              <span>Tick labels (pt)</span>
              <input type="number" min={3} max={24} step={0.5} value={tickPt(doc)}
                onChange={e => onChange({ ...doc, font: { ...doc.font, tickPt: num(e.target.value, 3, 24, 7) } })} />
            </label>
            <label className={styles.field}>
              <span>Letters (pt)</span>
              <input type="number" min={6} max={48} value={doc.letters.sizePt}
                onChange={e => onChange({ ...doc, letters: { ...doc.letters, sizePt: num(e.target.value, 6, 48, 16) } })} />
            </label>
            <label className={styles.field}>
              <span>Letter case</span>
              <select value={doc.letters.case}
                onChange={e => onChange({ ...doc, letters: { ...doc.letters, case: e.target.value as FigureDoc['letters']['case'] } })}>
                <option value="lower">a, b, c</option>
                <option value="upper">A, B, C</option>
                <option value="none">None</option>
              </select>
            </label>
            <label className={styles.field}>
              <span>Group gap (mm)</span>
              <input type="number" min={0} max={40} value={doc.groupGapMm}
                onChange={e => onChange({ ...doc, groupGapMm: num(e.target.value, 0, 40, 4) })} />
            </label>
          </div>

          <ColourEditor doc={doc} onChange={onChange}
            names={mergeEntries(Object.values(paneEntries)).map(e => e.name)} />
          <LabelEditor doc={doc} onChange={onChange} />

          <div className={styles.section}>
            <div className={styles.sectionTitle}>
              Groups
              <span className={styles.spacer} />
              <button className="btn btn-sm" onClick={() => {
                onChange({ ...doc, groups: [...doc.groups, newGroup(1, 2)] })
                setSel({ group: doc.groups.length, pane: null })
              }}>Add group</button>
            </div>
            {doc.groups.map((g, i) => (
              <div key={g.id} className={`${styles.groupRow} ${i === sel.group ? styles.selected : ''}`}
                onClick={() => setSel({ group: i, pane: null })}>
                <strong style={{ minWidth: 18 }}>{groupName(i)}</strong>
                <span>{g.rows} × {g.cols}</span>
                <span className={styles.spacer} />
                <button className="btn btn-sm" aria-label="Move group up" disabled={i === 0}
                  onClick={e => { e.stopPropagation(); moveGroup(i, -1) }}>↑</button>
                <button className="btn btn-sm" aria-label="Move group down" disabled={i === doc.groups.length - 1}
                  onClick={e => { e.stopPropagation(); moveGroup(i, 1) }}>↓</button>
                <button className="btn btn-sm btn-danger" aria-label="Remove group" disabled={doc.groups.length === 1}
                  onClick={e => {
                    e.stopPropagation()
                    if (!window.confirm(`Remove group ${groupName(i)} and its panes?`)) return
                    onChange({ ...doc, groups: doc.groups.filter((_, k) => k !== i) })
                    setSel({ group: Math.max(0, i - 1), pane: null })
                  }}>✕</button>
              </div>
            ))}
          </div>

          {group && (
            <div className={styles.section}>
              <div className={styles.sectionTitle}>Group {groupName(sel.group)}</div>
              <label className={styles.field}>
                <span>Title</span>
                <input type="text" value={group.title ?? ''} placeholder="None"
                  onChange={e => setGroup(sel.group, { ...group, title: e.target.value || null })} />
              </label>
              <label className={styles.field}>
                <span>Letter</span>
                <input type="text" value={group.label ?? ''} placeholder="Automatic" maxLength={4}
                  onChange={e => setGroup(sel.group, { ...group, label: e.target.value === '' ? null : e.target.value })} />
              </label>
              <div className={styles.field}>
                <span>Grid</span>
                <span className={styles.inline}>
                  <input type="number" min={1} max={6} value={group.rows} aria-label="Rows"
                    onChange={e => setGroup(sel.group, resizeGroup(group, num(e.target.value, 1, 6, 1), group.cols))} />
                  rows ×
                  <input type="number" min={1} max={6} value={group.cols} aria-label="Columns"
                    onChange={e => setGroup(sel.group, resizeGroup(group, group.rows, num(e.target.value, 1, 6, 1)))} />
                </span>
              </div>
              {group.cols > 1 && (
                <div className={styles.field}>
                  <span>Column widths</span>
                  <span className={styles.inline}>
                    {Array.from({ length: group.cols }, (_, c) => (
                      <input key={c} type="number" min={0.25} max={12} step={0.25} style={{ width: 52 }}
                        value={group.colWeights?.[c] ?? 1} aria-label={`Width of column ${c + 1}`}
                        onChange={e => {
                          const w = Array.from({ length: group.cols }, (_, k) => group.colWeights?.[k] ?? 1)
                          w[c] = num(e.target.value, 0.25, 12, 1)
                          setGroup(sel.group, { ...group, colWeights: w.every(x => x === 1) ? undefined : w })
                        }} />
                    ))}
                  </span>
                </div>
              )}
              {group.rows > 1 && (
                <div className={styles.field}>
                  <span>Row heights</span>
                  <span className={styles.inline}>
                    {Array.from({ length: group.rows }, (_, r) => (
                      <input key={r} type="number" min={0.25} max={12} step={0.05} style={{ width: 52 }}
                        value={group.rowWeights?.[r] ?? 1} aria-label={`Height of row ${r + 1}`}
                        onChange={e => {
                          const w = Array.from({ length: group.rows }, (_, k) => group.rowWeights?.[k] ?? 1)
                          w[r] = num(e.target.value, 0.25, 12, 1)
                          setGroup(sel.group, { ...group, rowWeights: w.every(x => x === 1) ? undefined : w })
                        }} />
                    ))}
                  </span>
                </div>
              )}
              <label className={styles.field}>
                <span>Height share</span>
                <input type="number" min={0.25} max={12} step={0.25} value={group.weight}
                  onChange={e => setGroup(sel.group, { ...group, weight: num(e.target.value, 0.25, 12, 1) })} />
              </label>
              <label className={styles.field}>
                <span>Legend</span>
                <select value={legendPlacement(group)}
                  onChange={e => setGroup(sel.group, { ...group, legend: e.target.value as LegendPlacement })}>
                  <option value="right">One, at the right</option>
                  <option value="top">One row, above the charts</option>
                  <option value="bottom">One row, below the charts</option>
                  <option value="each">On every pane</option>
                  <option value="none">None</option>
                </select>
              </label>
              {(legendPlacement(group) === 'top' || legendPlacement(group) === 'bottom') && (
                <label className={styles.field}>
                  <span>Italic legend</span>
                  <input type="checkbox" checked={group.legendItalic ?? false}
                    onChange={e => setGroup(sel.group, { ...group, legendItalic: e.target.checked || undefined })} />
                </label>
              )}
              {(legendPlacement(group) === 'top' || legendPlacement(group) === 'bottom') && (
                <label className={styles.field}>
                  <span>Bold legend</span>
                  <input type="checkbox" checked={group.legendBold ?? false}
                    onChange={e => setGroup(sel.group, { ...group, legendBold: e.target.checked || undefined })} />
                </label>
              )}
              <div className={styles.field}>
                <span>Background</span>
                <span className={styles.inline}>
                  <input type="checkbox" checked={group.background != null} aria-label="Tint the group background"
                    onChange={e => setGroup(sel.group, { ...group, background: e.target.checked ? '#e5e5e5' : null })} />
                  {group.background != null && (
                    <input type="color" value={group.background} aria-label="Background colour"
                      onChange={e => setGroup(sel.group, { ...group, background: e.target.value })} />
                  )}
                </span>
              </div>
            </div>
          )}

          {group && sel.pane != null && group.panes[sel.pane] && (
            <PaneEditor study={study} figureId={doc.id} runs={runs} source={source} pane={group.panes[sel.pane]}
              onChange={p => setPane(sel.group, sel.pane!, p)} />
          )}
          {group && sel.pane == null && (
            <p style={{ color: 'var(--color-muted-fg)', margin: 0 }}>Select a pane on the page to choose what it shows.</p>
          )}
        </div>

        <div ref={stageRef} className={styles.stage}>
          <div style={{ width: pageW * scale, height: pageH * scale }}>
            <div className={styles.page} style={{ width: pageW, height: pageH, transform: `scale(${scale})` }}>
              {doc.groups.map((g, gi) => {
                const gg = geo.groups[gi]
                const b = gg.box
                return (
                  <div key={g.id}>
                    {g.background && (
                      <div className={styles.groupBox} style={{
                        left: b.x * PX_PER_MM, top: b.y * PX_PER_MM, width: b.w * PX_PER_MM, height: b.h * PX_PER_MM,
                        background: g.background,
                      }} />
                    )}
                    {g.panes.map((p, pi) => {
                      const box = gg.panes[pi]
                      const w = box.w * PX_PER_MM, h = box.h * PX_PER_MM
                      const selected = sel.group === gi && sel.pane === pi
                      return (
                        <div key={pi} role="button" tabIndex={0}
                          aria-label={`Group ${groupName(gi)}, pane ${pi + 1}`}
                          className={`${styles.pane} ${selected ? styles.selected : ''}`}
                          style={{ left: box.x * PX_PER_MM, top: box.y * PX_PER_MM, width: w, height: h }}
                          onClick={() => setSel({ group: gi, pane: pi })}
                          onKeyDown={e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); setSel({ group: gi, pane: pi }) } }}>
                          <PaneContent study={study} figureId={doc.id} pane={p} doc={doc}
                            showLegend={showsLegend(g, pi)} w={w} h={h} version={version}
                            onEntries={e => reportEntries(`${g.id}:${pi}`, e)} />
                        </div>
                      )
                    })}
                    {gg.title && (
                      <div style={{ position: 'absolute', left: gg.title.x * PX_PER_MM, top: gg.title.y * PX_PER_MM,
                                    width: gg.title.w * PX_PER_MM, height: gg.title.h * PX_PER_MM,
                                    lineHeight: `${gg.title.h * PX_PER_MM}px`, textAlign: 'center', color: '#444',
                                    fontFamily: fontFamily(doc), fontSize: titlePt(doc) * PX_PER_PT, pointerEvents: 'none',
                                    fontWeight: doc.font.titleBold ? 'bold' : undefined }}>
                        {g.title}
                      </div>
                    )}
                    {gg.legend && (
                      <GroupLegend doc={doc} box={gg.legend} italic={g.legendItalic} bold={g.legendBold}
                        entries={mergeEntries(g.panes.map((p, pi) => p.chart ? paneEntries[`${g.id}:${pi}`] ?? [] : []))} />
                    )}
                    {gg.letter && <span className={styles.letter} style={{
                      left: b.x * PX_PER_MM + doc.letters.sizePt * PX_PER_PT * 0.15,
                      top: b.y * PX_PER_MM,
                      fontSize: doc.letters.sizePt * PX_PER_PT,
                      fontFamily: doc.font.family === 'Arial' ? 'Arial, Helvetica, sans-serif' : '"Times New Roman", Times, serif',
                    }}>{gg.letter}</span>}
                  </div>
                )
              })}
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}

/** Colours by legend name, applied to every pane of the figure. */
function ColourEditor({ doc, names, onChange }: { doc: FigureDoc; names: string[]; onChange: (doc: FigureDoc) => void }) {
  const colours = doc.colours ?? {}
  const [adding, setAdding] = useState('')
  const set = (next: Record<string, string>) => onChange({ ...doc, colours: next })
  const unused = names.filter(n => !(n in colours))
  return (
    <div className={styles.section}>
      <div className={styles.sectionTitle}>Colours</div>
      {Object.entries(colours).map(([name, c]) => (
        <div key={name} className={styles.field}>
          <span title={name} style={{ overflow: 'hidden', textOverflow: 'ellipsis' }}>{name}</span>
          <span className={styles.inline}>
            <input type="color" value={toHex(c)} aria-label={`Colour for ${name}`}
              onChange={e => set({ ...colours, [name]: e.target.value })} />
            <button className="btn btn-sm" aria-label={`Use the chart's colour for ${name}`}
              onClick={() => { const { [name]: _, ...rest } = colours; set(rest) }}>✕</button>
          </span>
        </div>
      ))}
      {unused.length > 0 && (
        <div className={styles.inline}>
          <select value={adding} onChange={e => setAdding(e.target.value)} aria-label="Series to colour">
            <option value="">Choose a series…</option>
            {unused.map(n => <option key={n} value={n}>{n}</option>)}
          </select>
          <button className="btn btn-sm" disabled={!adding}
            onClick={() => { set({ ...colours, [adding]: '#888888' }); setAdding('') }}>Set colour</button>
        </div>
      )}
    </div>
  )
}

/** Display text for run, group and category labels, applied to every pane. */
function LabelEditor({ doc, onChange }: { doc: FigureDoc; onChange: (doc: FigureDoc) => void }) {
  const labels = doc.labels ?? {}
  const [from, setFrom] = useState('')
  const [to, setTo] = useState('')
  const set = (next: Record<string, string>) => onChange({ ...doc, labels: next })
  return (
    <div className={styles.section}>
      <div className={styles.sectionTitle}>Labels</div>
      {Object.entries(labels).map(([name, shown]) => (
        <div key={name} className={styles.field}>
          <span title={name} style={{ overflow: 'hidden', textOverflow: 'ellipsis' }}>{name}</span>
          <span className={styles.inline} style={{ flexWrap: 'nowrap' }}>
            <input type="text" value={shown} aria-label={`Shown text for ${name}`} style={{ flex: 1, minWidth: 0 }}
              onChange={e => set({ ...labels, [name]: e.target.value })} />
            <button className="btn btn-sm" aria-label={`Show ${name} as it is`}
              onClick={() => { const { [name]: _, ...rest } = labels; set(rest) }}>✕</button>
          </span>
        </div>
      ))}
      <div className={styles.inline} style={{ flexWrap: 'nowrap' }}>
        <input type="text" value={from} placeholder="Label, e.g. Large Intestine" aria-label="Label to replace"
          style={{ flex: 1, minWidth: 0, font: 'inherit' }} onChange={e => setFrom(e.target.value)} />
        <input type="text" value={to} placeholder="Shown as" aria-label="Text to show instead"
          style={{ flex: 1, minWidth: 0, font: 'inherit' }} onChange={e => setTo(e.target.value)} />
        <button className="btn btn-sm" disabled={!from.trim()}
          onClick={() => { set({ ...labels, [from.trim()]: to }); setFrom(''); setTo('') }}>Add</button>
      </div>
      <span style={{ color: 'var(--color-muted-fg)' }}>Underscores already show as spaces; name the label with spaces.</span>
    </div>
  )
}

function toHex(c: string): string {
  const p = parseColour(c)
  return p ? '#' + p.slice(0, 3).map(v => Math.round(v).toString(16).padStart(2, '0')).join('') : '#888888'
}
