// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { api } from '../api/client'
import { extractAlphaPanel } from '../components/alphaMetrics'
import { textWidth } from './legend'
import { PX_PER_PT, tickPt, titlePt, type AlphaMetricKey, type ChartSpec, type FigureDoc, type FigurePane } from './types'

export interface PlotSpec {
  data: Record<string, unknown>[]
  layout: Record<string, unknown>
}

// One request per distinct body, shared by every pane that needs it.
const requests = new Map<string, Promise<unknown>>()

export function clearChartCache() { requests.clear() }

function once(key: string, make: () => Promise<unknown>): Promise<unknown> {
  let p = requests.get(key)
  if (!p) {
    p = make()
    p.catch(() => requests.delete(key))
    requests.set(key, p)
  }
  return p
}

const legendKey = (t: Record<string, unknown>) => String(t.legendgroup ?? t.name ?? '')

// The alpha figure lists its legend on the first metric only, so a pane showing
// another metric takes over those entries.
function alphaPanel(fig: unknown, metric: AlphaMetricKey): unknown {
  const all = ((fig as { data?: Record<string, unknown>[] })?.data ?? [])
  const listed = new Set(all.filter(t => t.showlegend !== false && legendKey(t)).map(legendKey))
  const panel = extractAlphaPanel(fig, metric) as { data?: Record<string, unknown>[] }
  const seen = new Set<string>()
  const data = (panel.data ?? []).map(t => {
    const k = legendKey(t)
    const show = listed.has(k) && !seen.has(k)
    if (show) seen.add(k)
    return { ...t, showlegend: show }
  })
  return { ...panel, data }
}

/** The chart as the backend draws it, before any pane styling. */
export async function fetchChart(study: string, spec: ChartSpec): Promise<unknown> {
  const runs = spec.runs.map(r => ({ ...r, source: spec.source }))
  if (spec.kind === 'alpha' || spec.kind === 'nmds') {
    const body = { table: spec.table, runs, colFilters: {}, aggregate: spec.aggregate }
    const key = `${spec.kind}|${study}|${JSON.stringify(body)}`
    const fig = await once(key, () => spec.kind === 'alpha'
      ? api.analysis.compareAlpha(study, body)
      : api.analysis.nmds(study, body))
    return spec.kind === 'alpha' ? alphaPanel(fig, spec.metric) : fig
  }
  const body = {
    table: spec.table, tag: spec.tag, value: spec.value, relative: spec.relative, mode: spec.mode,
    keep_empty: spec.keepEmpty, ...(spec.tag === 'rank' ? { top_n: spec.topN } : {}),
  }
  if (spec.facet) {
    const req = { ...body, runs, rows: spec.facet.rows, cols: spec.facet.cols }
    return once(`facet|${study}|${JSON.stringify(req)}`, () => api.analysis.chartFacet(study, req))
  }
  const req = { ...body, runs, subgroup: null }
  return once(`chart|${study}|${JSON.stringify(req)}`, () => api.analysis.chartCompare(study, req))
}

const axisKeys = (layout: Record<string, unknown>, letter: 'x' | 'y') =>
  Object.keys(layout).filter(k => new RegExp(`^${letter}axis\\d*$`).test(k))

const titleText = (t: unknown): string =>
  typeof t === 'string' ? t : (t as { text?: string } | undefined)?.text ?? ''

/** [r, g, b, a] of a hex or rgb()/rgba() colour, or null. */
export function parseColour(c: string): [number, number, number, number] | null {
  const hex = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(c.trim())
  if (hex) {
    const h = hex[1].length === 3 ? hex[1].split('').map(x => x + x).join('') : hex[1]
    return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16), 1]
  }
  const rgb = /^rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*(?:,\s*([\d.]+)\s*)?\)$/i.exec(c.trim())
  return rgb ? [Number(rgb[1]), Number(rgb[2]), Number(rgb[3]), rgb[4] == null ? 1 : Number(rgb[4])] : null
}

const rgba = ([r, g, b]: number[], a: number) => `rgba(${Math.round(r)},${Math.round(g)},${Math.round(b)},${a})`

/** Recolour the traces named in `colours`: box fills keep the chart's transparency,
 *  and box outlines and points take a darker shade. */
function recolour(t: Record<string, unknown>, colours: Record<string, string>): Record<string, unknown> {
  const c = typeof t.name === 'string' ? parseColour(colours[t.name] ?? '') : null
  if (!c) return t
  const dark = rgba(c.map(v => v * 0.7), 1)
  const marker = { ...(t.marker as object | undefined) } as Record<string, unknown>
  if (t.type === 'box') {
    const alpha = parseColour(String(t.fillcolor ?? ''))?.[3] ?? 0.8
    return { ...t, fillcolor: rgba(c, alpha), line: { ...(t.line as object | undefined), color: dark },
             marker: { ...marker, color: dark } }
  }
  const outline = { ...(marker.line as object | undefined), color: dark }
  return { ...t, marker: { ...marker, color: rgba(c, c[3]), ...(t.type === 'scatter' ? { line: outline } : {}) } }
}

/** The chart sized, typeset and titled for its pane. */
export function paneFigure(raw: unknown, pane: FigurePane, doc: FigureDoc, opts: {
  showLegend: boolean
  widthPx: number
  heightPx: number
}): PlotSpec {
  const fig = JSON.parse(JSON.stringify(raw)) as PlotSpec
  const base = doc.font.sizePt * PX_PER_PT
  const tick = tickPt(doc) * PX_PER_PT
  const style = pane.style ?? {}
  const family = doc.font.family === 'Arial' ? 'Arial, Helvetica, sans-serif' : '"Times New Roman", Times, serif'
  const L = fig.layout ?? {}

  const title = pane.title ?? titleText(L.title)
  const hasTitle = title.trim().length > 0
  // A title wider than its pane shrinks to fit, clear of the group letter.
  const titleRoom = opts.widthPx - 2 * doc.letters.sizePt * PX_PER_PT
  const titleSize = Math.min(titlePt(doc) * PX_PER_PT,
    titlePt(doc) * PX_PER_PT * titleRoom / Math.max(1, textWidth(title.replace(/<[^>]+>/g, ''), titlePt(doc) * PX_PER_PT, family)))

  for (const k of axisKeys(L, 'x').concat(axisKeys(L, 'y'))) {
    const ax = { ...(L[k] as Record<string, unknown>) }
    ax.automargin = true
    ax.tickfont = { ...(ax.tickfont as object | undefined), size: tick }
    if (k.startsWith('x') && style.tickAngle != null) ax.tickangle = style.tickAngle
    if (k.startsWith('x') && style.allTicks) { ax.tickmode = 'linear'; ax.dtick = 1 }
    if (ax.title) ax.title = { ...(typeof ax.title === 'string' ? { text: ax.title } : ax.title as object), font: { size: base } }
    L[k] = ax
  }
  if (pane.xTitle != null) {
    const xs = axisKeys(L, 'x')
    const titled = xs.filter(k => titleText((L[k] as Record<string, unknown>).title))
    for (const k of titled.length ? titled : xs.slice(0, 1))
      L[k] = { ...(L[k] as object), title: { text: pane.xTitle, font: { size: base } } }
  }
  // Grid charts head each row and column with "run: X"; a figure shows just X.
  const facetLabel = /(^|>)(run|group|subgroup): /
  const isGrid = L.grid != null ||
    (Array.isArray(L.annotations) && (L.annotations as { text?: string }[]).some(a => facetLabel.test(a.text ?? '')))
  if (Array.isArray(L.annotations)) {
    L.annotations = (L.annotations as Record<string, unknown>[]).map(a => ({
      ...a,
      ...(typeof a.text === 'string' && facetLabel.test(a.text)
        ? { text: a.text.replace(facetLabel, '$1').replace(/_/g, ' ') } : {}),
      font: { ...(a.font as object | undefined), size: base * 0.9 },
    }))
  }

  if (L.grid && style.gridGap != null) L.grid = { ...(L.grid as object), xgap: style.gridGap }
  if (L.grid && style.gridRowGap != null) L.grid = { ...(L.grid as object), ygap: style.gridRowGap }

  fig.data = (fig.data ?? []).map(raw => {
    let t = recolour(raw.type === 'box' && raw.width == null ? { ...raw, width: 0.9 } : raw, doc.colours ?? {})
    const lineOnly = t.type === 'scatter' && t.mode === 'lines'
    if (lineOnly && style.lineOpacity != null) {
      const c = parseColour(String((t.line as Record<string, unknown> | undefined)?.color ?? '')) ?? [70, 70, 70, 1]
      t = { ...t, line: { ...(t.line as object | undefined), color: rgba(c, style.lineOpacity) } }
    }
    const hasPoints = t.type === 'box' || (t.type === 'scatter' && String(t.mode ?? 'markers').includes('markers'))
    if (hasPoints && style.pointBorder != null) {
      const m = (t.marker ?? {}) as Record<string, unknown>
      t = { ...t, marker: { ...m, line: { ...(m.line as object | undefined), width: style.pointBorder } } }
    }
    if (hasPoints && style.pointOpacity != null) t = { ...t, marker: { ...(t.marker as object | undefined), opacity: style.pointOpacity } }
    return t
  })
  fig.layout = {
    ...L,
    width: opts.widthPx, height: opts.heightPx, autosize: false,
    font: { ...(L.font as object | undefined), family, size: base },
    title: hasTitle ? { text: title, font: { size: titleSize }, x: 0.5, xanchor: 'center', y: 1, yanchor: 'top', yref: 'container', pad: { t: 4 } } : { text: '' },
    showlegend: opts.showLegend,
    legend: { ...(L.legend as object | undefined), font: { size: base } },
    margin: { l: 8, r: 8, t: (hasTitle ? titleSize * 2.1 : 8) + (isGrid ? base * 1.6 : 0), b: 8, pad: 0 },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
  }
  return fig
}
