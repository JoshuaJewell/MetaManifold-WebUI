// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

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

/** Gap between horizontal bars, as a share of each category's height. */
const HORIZONTAL_BARGAP = 0.08

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
  // A title wider than its pane shrinks to fit, clear of the group letter when letters are on.
  const titleRoom = opts.widthPx - (doc.letters.case === 'none' ? 0 : 2 * doc.letters.sizePt * PX_PER_PT)
  const titleSize = Math.min(titlePt(doc) * PX_PER_PT,
    titlePt(doc) * PX_PER_PT * titleRoom / Math.max(1, textWidth(title.replace(/<[^>]+>/g, ''), titlePt(doc) * PX_PER_PT, family, false, doc.font.titleBold)))

  // Runs of grouped studies are labelled "group/run". A group every label shares is
  // dropped; otherwise pooled bars get a two-level axis, run under group. Underscores show as
  // spaces, and the figure's labels replace any label they name.
  const shown = (label: string) => doc.labels?.[label] ?? label
  const isPair = (v: unknown) => typeof v === 'string' && v.split('/').length === 2
  const allPairs = (ts: Record<string, unknown>[]) => ts.length > 0 &&
    ts.every(t => Array.isArray(t.x) && t.x.length > 0 && (t.x as unknown[]).every(isPair))
  const bars = (fig.data ?? []).filter(t => t.type === 'bar')
  const others = (fig.data ?? []).filter(t => t.type === 'box' || t.type === 'scatter')
  const groupsOf = (ts: Record<string, unknown>[]) => new Set(ts.flatMap(t => (t.x as string[]).map(v => v.split('/')[0])))
  if (allPairs(bars) && groupsOf(bars).size === 1) {
    fig.data = (fig.data ?? []).map(t => t.type !== 'bar' ? t
      : { ...t, x: (t.x as string[]).map(v => shown(v.split('/')[1].replace(/_/g, ' '))) })
  } else if (allPairs(bars)) {
    fig.data = (fig.data ?? []).map(t => {
      if (t.type !== 'bar') return t
      const parts = (t.x as string[]).map(v => v.split('/').map(p => shown(p.replace(/_/g, ' '))))
      return { ...t, x: [parts.map(p => p[0]), parts.map(p => p[1])] }
    })
    for (const k of axisKeys(L, 'x')) L[k] = { ...(L[k] as object), type: 'multicategory', tickangle: 0 }
  } else if (!bars.length && allPairs(others)) {
    const groups = groupsOf(others)
    const categories = new Set(others.flatMap(t => t.x as string[]))
    // A label wider than its slot along the axis breaks at its spaces.
    const slot = opts.widthPx * 0.8 / Math.max(1, categories.size)
    const label = (v: string) => {
      const text = shown((groups.size === 1 ? v.split('/')[1] : v).replace(/_/g, ' '))
      return textWidth(text, tickPt(doc) * PX_PER_PT, family) > slot ? text.replace(/ /g, '<br>') : text
    }
    fig.data = (fig.data ?? []).map(t => others.includes(t) ? { ...t, x: (t.x as string[]).map(label) } : t)
  }

  // Segment labels: each series' name inside its segments; ones too small to read are hidden.
  // Genus and species names are set in italics.
  if (style.labelSegments && bars.length) {
    const italic = pane.chart?.kind === 'composition' && pane.chart.tag === 'rank' && /^(genus|species)$/i.test(pane.chart.value)
    fig.data = (fig.data ?? []).map(t => {
      if (t.type !== 'bar') return t
      const name = String(t.name ?? '').replace(/_/g, ' ')
      const n = (Array.isArray(t.y) ? t.y : []).length
      return { ...t, text: Array(n).fill(italic && name !== 'Other' ? `<i>${name}</i>` : name), texttemplate: '%{text}',
               textposition: 'inside', insidetextanchor: 'middle', constraintext: 'inside', textangle: 'auto' }
    })
    L.uniformtext = { mode: 'hide', minsize: Math.min(tickPt(doc), 6) * PX_PER_PT }
  }

  // Horizontal: values run along x and categories down y, first at the top. Bars and
  // boxes turn, and the lines joining paired samples turn with them.
  if (style.horizontal && (bars.length || others.length)) {
    fig.data = (fig.data ?? []).map(t => t.type !== 'bar' && t.type !== 'box' && t.type !== 'scatter' ? t
      : { ...t, x: t.y, y: t.x, ...(t.type === 'scatter' ? {} : { orientation: 'h' }),
          xaxis: t.yaxis?.toString().replace('y', 'x'), yaxis: t.xaxis?.toString().replace('x', 'y') })
    const xs = axisKeys(L, 'x'), ys = axisKeys(L, 'y')
    const old = Object.fromEntries([...xs, ...ys].map(k => [k, L[k]]))
    for (const k of xs) delete L[k]
    for (const k of ys) delete L[k]
    for (const k of xs) L[k.replace('x', 'y')] = { ...(old[k] as object), autorange: 'reversed', tickangle: 0 }
    for (const k of ys) L[k.replace('y', 'x')] = { ...(old[k] as object), tickangle: 0 }
    // Thick bars leave room for labels inside them.
    if (bars.length) L.bargap = HORIZONTAL_BARGAP
  }

  // Plotly's automatic margin leaves no room for a two-level category axis down the
  // side, so the labels would spill out of the pane: measure them instead.
  // Horizontal charts measure their category labels too, so the plot's width is known
  // and category bands can reach the pane's edges.
  let left = 8
  const tickPx = tickPt(doc) * PX_PER_PT
  const widest = (labels: unknown[] = []) => Math.max(0, ...labels.flatMap(l =>
    String(l).replace(/<(?!br)[^>]+>/g, '').split('<br>').map(line => textWidth(line, tickPx, family))))
  const yMulti = (L.yaxis as Record<string, unknown> | undefined)?.type === 'multicategory'
  const turned = (fig.data ?? []).find(t => t.orientation === 'h' && Array.isArray(t.y))
  if (yMulti && !style.hideYTicks) {
    const levels = (turned?.y ?? []) as string[][]
    left = widest(levels[0]) + widest(levels[1]) + tickPx * 2.5
  } else if (turned && style.horizontal && !style.hideYTicks) {
    left = widest(turned.y as unknown[]) + tickPx * 1.2
    for (const k of axisKeys(L, 'y')) L[k] = { ...(L[k] as object), automargin: false }
  }

  for (const k of axisKeys(L, 'x').concat(axisKeys(L, 'y'))) {
    const ax = { ...(L[k] as Record<string, unknown>) }
    ax.automargin = true
    ax.tickfont = { ...(ax.tickfont as object | undefined), size: tick }
    if (k.startsWith('x') && style.tickAngle != null) ax.tickangle = style.tickAngle
    if (k.startsWith('x') && style.allTicks) { ax.tickmode = 'linear'; ax.dtick = 1 }
    if (k.startsWith('x') ? style.hideXTicks : style.hideYTicks) ax.showticklabels = false
    if (ax.title) ax.title = { ...(typeof ax.title === 'string' ? { text: ax.title } : ax.title as object), font: { size: base } }
    L[k] = ax
  }
  if (pane.xTitle != null) {
    const xs = axisKeys(L, 'x')
    const titled = xs.filter(k => titleText((L[k] as Record<string, unknown>).title))
    for (const k of titled.length ? titled : xs.slice(0, 1))
      L[k] = { ...(L[k] as object), title: { text: pane.xTitle, font: { size: base } } }
  }
  if (pane.yTitle != null) {
    for (const k of axisKeys(L, 'y'))
      if (titleText((L[k] as Record<string, unknown>).title))
        L[k] = { ...(L[k] as object), title: { text: pane.yTitle, font: { size: base } } }
  }
  // Grid charts head each row and column with "run: X"; a figure shows just X.
  const facetLabel = /(^|>)(run|group|subgroup): /
  const gridCells = L.grid ? ((L.grid as { rows?: number }).rows ?? 1) * ((L.grid as { columns?: number }).columns ?? 1) : 0
  const isGrid = gridCells > 1 ||
    (Array.isArray(L.annotations) && (L.annotations as { text?: string }[]).some(a => facetLabel.test(a.text ?? '')))
  if (style.hideNotes && Array.isArray(L.annotations))
    L.annotations = (L.annotations as { text?: string }[]).filter(a => facetLabel.test(a.text ?? ''))
  // A header wider than its panel breaks at its spaces.
  const grid = L.grid as { rows?: number; columns?: number } | undefined
  let wrappedHeader = false
  if (Array.isArray(L.annotations)) {
    L.annotations = (L.annotations as Record<string, unknown>[])
      .filter(a => !(style.hideRowHeaders && a.textangle === 90 && typeof a.text === 'string' && facetLabel.test(a.text)))
      .map(a => {
      if (typeof a.text !== 'string' || !facetLabel.test(a.text))
        return { ...a, font: { ...(a.font as object | undefined), size: base * 0.9 } }
      let text = a.text.replace(facetLabel, '$1').replace(/_/g, ' ')
        .replace(/(<b>)?([^<]+)(<\/b>)?/, (m, o, label, c) => doc.labels?.[label] != null ? `${o ?? ''}${doc.labels[label]}${c ?? ''}` : m)
      const vertical = a.textangle === 90
      const room = vertical ? opts.heightPx / (grid?.rows ?? 1) : opts.widthPx / (grid?.columns ?? 1)
      const plain = text.replace(/<[^>]+>/g, '')
      // Header text is bold, so it is measured bold.
      const width = (t: string) => textWidth(t, base * 0.9, family, false, true)
      if (plain.includes(' ') && width(plain) > room * 0.9) {
        text = text.replace(/ /g, '<br>')
        if (!vertical) wrappedHeader = true
      }
      // A line still too wide shrinks, down to the tick label size.
      const widest = Math.max(...text.replace(/<b>|<\/b>/g, '').split('<br>').map(width))
      const size = Math.max(tick, base * 0.9 * Math.min(1, room * 0.9 / widest))
      return { ...a, text, font: { ...(a.font as object | undefined), size } }
    })
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
    if (t.type === 'box' && style.fillOpacity != null) {
      const c = parseColour(String(t.fillcolor ?? ''))
      if (c) t = { ...t, fillcolor: rgba(c, style.fillOpacity) }
    }
    return t
  })
  // The value axis runs along x once turned, along y otherwise. A fixed range shows
  // three ticks: its ends and its middle.
  const valueAxis = style.horizontal ? 'xaxis' : 'yaxis'
  if (style.valueRange) {
    const [a, b] = style.valueRange
    const step = style.valueStep && style.valueStep > 0 ? style.valueStep : (b - a) / 2
    const ticks: number[] = []
    for (let v = a; v <= b + step * 1e-6 && ticks.length < 50; v += step) ticks.push(Math.round(v * 1e9) / 1e9)
    L[valueAxis] = { ...(L[valueAxis] as object), range: [a, b], autorange: false, tickmode: 'array', tickvals: ticks }
  }
  if (style.grid === 'none') L[valueAxis] = { ...(L[valueAxis] as object), showgrid: false }
  if (style.grid === 'strong') L[valueAxis] = { ...(L[valueAxis] as object), showgrid: true, gridcolor: '#a8a8a8', gridwidth: 1 }

  // Categories in order of first appearance, as Plotly places them.
  const catAxis = style.horizontal ? 'y' : 'x'
  const cats = [...new Set((fig.data ?? []).filter(t => t.type !== 'scatter')
    .flatMap(t => Array.isArray(t[catAxis]) ? (t[catAxis] as unknown[]).map(String) : []))]

  // A comparison bracket beside the first two categories, labelled, for a test the
  // chart could not run itself or whose result is given by hand.
  let right = 8
  const compare = style.compareLabel?.trim()
  const bracketGap = 3, bracketTick = 3
  if (compare && cats.length >= 2) {
    const labelPx = textWidth(compare, tickPx, family, false, true)
    if (style.horizontal) right = bracketGap + bracketTick + 3 + labelPx + 2
  }
  // Horizontal charts set every margin themselves: Plotly's automatic margins move each
  // plot by a pixel or two with its tick labels, so panes side by side would not line up.
  const xTicksShown = !(style.horizontal ? style.hideXTicks : style.hideYTicks)
  // Horizontal charts keep their margins tight: a band's overhang, and the tick labels.
  const bandPad = style.categoryBands ? 4 : 0
  const xTitled = titleText((L.xaxis as Record<string, unknown> | undefined)?.title).trim().length > 0
  const bottom = style.horizontal ? bandPad + 1 + (xTicksShown ? tickPx * 1.25 : 0) + (xTitled ? base * 1.15 : 0) : 8
  if (style.horizontal)
    for (const k of axisKeys(L, 'x').concat(axisKeys(L, 'y'))) {
      const ax = L[k] as Record<string, unknown>
      // The title sits just under the tick labels.
      const t = ax.title && typeof ax.title === 'object' ? { ...(ax.title as object), standoff: 1 } : ax.title
      L[k] = { ...ax, automargin: false, ...(t ? { title: t } : {}) }
    }
  const plotW = Math.max(1, opts.widthPx - left - right)
  const topPx = (hasTitle ? titlePt(doc) * PX_PER_PT * (style.horizontal ? 1.45 : 2.1) : style.horizontal ? bandPad + 1 : 8)
    + (isGrid ? base * (wrappedHeader ? 2.7 : 1.6) : 0)
  const plotH = Math.max(1, opts.heightPx - topPx - bottom)
  if (compare && cats.length >= 2) {
    const yref = style.horizontal ? (L.yaxis ? 'y' : 'y') : 'paper'
    const text = `<b>${compare.replace(/</g, '&lt;').replace(/>/g, '&gt;')}</b>`
    const line = { color: '#444444', width: 1 }
    const shapes = style.horizontal
      ? (() => {
          const x = 1 + bracketGap / plotW, t = bracketTick / plotW
          return [{ type: 'line', xref: 'paper', yref, x0: x, x1: x, y0: 0, y1: 1, line },
                  { type: 'line', xref: 'paper', yref, x0: x - t, x1: x, y0: 0, y1: 0, line },
                  { type: 'line', xref: 'paper', yref, x0: x - t, x1: x, y0: 1, y1: 1, line }]
        })()
      : (() => {
          const y = 1 + bracketGap / plotH, t = bracketTick / plotH
          return [{ type: 'line', xref: 'x', yref: 'paper', x0: 0, x1: 1, y0: y, y1: y, line },
                  { type: 'line', xref: 'x', yref: 'paper', x0: 0, x1: 0, y0: y - t, y1: y, line },
                  { type: 'line', xref: 'x', yref: 'paper', x0: 1, x1: 1, y0: y - t, y1: y, line }]
        })()
    const note = style.horizontal
      ? { text, xref: 'paper', yref, x: 1 + (bracketGap + 3) / plotW, y: 0.5, xanchor: 'left', yanchor: 'middle', showarrow: false, font: { size: tickPx } }
      : { text, xref: 'x', yref: 'paper', x: 0.5, y: 1 + (bracketGap + 2) / plotH, xanchor: 'center', yanchor: 'bottom', showarrow: false, font: { size: tickPx } }
    L.shapes = [...((L.shapes as unknown[]) ?? []), ...shapes]
    L.annotations = [...((L.annotations as unknown[]) ?? []), note]
  }

  // Horizontal boxes leave the same gap between them as horizontal bars do. A box's outline
  // is stroked half outside it, so the box is narrowed by its line width.
  if (style.horizontal && cats.length) {
    const catPx = plotH / cats.length
    fig.data = (fig.data ?? []).map(t => {
      if (t.type !== 'box') return t
      const lw = Number((t.line as { width?: number } | undefined)?.width ?? 2)
      return { ...t, width: Math.max(0.1, 1 - HORIZONTAL_BARGAP - lw / catPx) }
    })
  }

  // Faint bands behind the categories named in the figure's colours, a little inside the
  // pane's edges and apart from each other, so neighbouring panes' bands read as rows.
  if (style.categoryBands) {
    const colourOf = (label: string) => {
      const plain = label.replace(/<br>/g, ' ').replace(/<[^>]+>/g, '').toLowerCase()
      const hit = Object.entries(doc.colours ?? {}).find(([k]) => k.replace(/_/g, ' ').toLowerCase() === plain)
      return hit ? parseColour(hit[1]) : null
    }
    // Horizontal bands run edge to edge across the pane, so they join the next pane's, meet
    // each other with no gap, and the outer ones overhang the plot by `bandPad` px.
    const pad = bandPad / (plotH / Math.max(1, cats.length))
    const shapes = cats.flatMap((c, i) => {
      const col = colourOf(c)
      if (!col) return []
      return [style.horizontal
        ? { type: 'rect', layer: 'below', xref: 'paper', yref: 'y', x0: -left / plotW, x1: 1 + right / plotW,
            y0: i - 0.5 - (i === 0 ? pad : 0), y1: i + 0.5 + (i === cats.length - 1 ? pad : 0),
            fillcolor: rgba(col, 0.28), line: { width: 0 } }
        : { type: 'rect', layer: 'below', xref: 'x', yref: 'paper', x0: i - 0.5, x1: i + 0.5,
            y0: 0, y1: 1, fillcolor: rgba(col, 0.28), line: { width: 0 } }]
    })
    // Bands go under everything else drawn on the chart.
    L.shapes = [...shapes, ...((L.shapes as unknown[]) ?? [])]
  }

  fig.layout = {
    ...L,
    width: opts.widthPx, height: opts.heightPx, autosize: false,
    font: { ...(L.font as object | undefined), family, size: base },
    title: hasTitle ? { text: doc.font.titleBold ? `<b>${title}</b>` : title, font: { size: titleSize }, x: 0.5, xref: 'paper', xanchor: 'center', y: 1, yanchor: 'top', yref: 'container', pad: { t: style.horizontal ? 1 : 4 } } : { text: '' },
    showlegend: opts.showLegend,
    legend: { ...(L.legend as object | undefined), font: { size: base } },
    // The margin is the full title size's even when a title shrinks, so neighbouring panes line up.
    margin: { l: left, r: right, t: topPx, b: bottom, pad: 0 },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
  }
  return fig
}
