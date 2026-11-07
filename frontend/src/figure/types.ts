// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import type { AnnotationSource, ComparisonRunSpec, FacetDimension } from '../api/types'

export type AlphaMetricKey = 'richness' | 'shannon' | 'simpson'

interface ChartBase {
  source: AnnotationSource
  runs: ComparisonRunSpec[]
  table: string
  aggregate: boolean
}

export type ChartSpec =
  | (ChartBase & { kind: 'alpha'; metric: AlphaMetricKey })
  | (ChartBase & { kind: 'nmds' })
  | (ChartBase & {
      kind: 'composition'
      tag: 'rank' | 'category'
      value: string
      topN: number
      relative: boolean
      mode: 'stacked' | 'grouped'
      keepEmpty: boolean
      /** A row × column grid of panels; null for one chart. */
      facet: { rows: FacetDimension; cols: FacetDimension } | null
    })

export type ChartKind = ChartSpec['kind']

/** A picture uploaded to the figure, drawn inside the pane. */
export interface PaneImage {
  /** Stored file name, from the figure's image upload. */
  file: string
  /** contain shows all of it; cover fills the pane and crops the overflow. */
  fit: 'contain' | 'cover'
}

/** One line of a text pane. *word* sets italics and **word** bold. */
export interface TextLine {
  text: string
  /** Size in pt; absent uses the figure's text size. */
  sizePt?: number
  colour?: string
}

export interface PaneText {
  lines: TextLine[]
  align: 'left' | 'center' | 'right'
  valign: 'top' | 'middle' | 'bottom'
  /** Turned a quarter: 90 reads top to bottom, -90 bottom to top. */
  rotate?: 0 | 90 | -90
  /** Moves the text right and down by this many mm, to line it up with something beside it. */
  offsetMm?: [number, number]
}

/** The box a text is laid out in, and the SVG transform that turns it into a w × h pane. */
export function textFrame(t: PaneText, w: number, h: number): { w: number; h: number; transform: string | undefined } {
  const [dx, dy] = (t.offsetMm ?? [0, 0]).map(v => v * PX_PER_MM)
  const shift = dx || dy ? `translate(${dx},${dy}) ` : ''
  if (t.rotate === -90) return { w: h, h: w, transform: `${shift}translate(0,${h}) rotate(-90)` }
  if (t.rotate === 90) return { w: h, h: w, transform: `${shift}translate(${w},0) rotate(90)` }
  return { w, h, transform: shift ? shift.trim() : undefined }
}

export interface FigurePane {
  chart: ChartSpec | null
  /** Shown when the pane has no chart. */
  image?: PaneImage | null
  /** Shown when the pane has no chart or image; under a picture, its caption. */
  text?: PaneText | null
  /** Replaces the chart's own title; empty hides it, null keeps it. */
  title: string | null
  /** Text for the x-axis title; null keeps the chart's own. */
  xTitle: string | null
  /** Text for every y-axis title; empty hides them, absent keeps the chart's own. */
  yTitle?: string | null
  style?: PaneStyle
}

/** Adjustments to the chart's own look; an absent key keeps the chart's. */
export interface PaneStyle {
  /** x tick label angle in degrees. */
  tickAngle?: number
  /** Label every x category instead of letting Plotly skip some. */
  allTicks?: boolean
  /** Leave out the x or y tick labels, as for a grid of many samples. */
  hideXTicks?: boolean
  hideYTicks?: boolean
  /** Opacity of line-only traces, such as the lines joining paired samples. */
  lineOpacity?: number
  /** Outline width of sample points, in px. */
  pointBorder?: number
  /** Opacity of box fills. */
  fillOpacity?: number
  /** Opacity of sample points. */
  pointOpacity?: number
  /** Draw bars or boxes horizontally, categories down the side. */
  horizontal?: boolean
  /** Write each series' name inside its bar segments, where it fits, in place of a legend. */
  labelSegments?: boolean
  /** Leave out the chart's own notes, such as test results or why a test did not run. */
  hideNotes?: boolean
  /** Leave this pane's series out of the group's legend. */
  hideFromLegend?: boolean
  /** Gridlines along the value axis: none, or darker so they show over category bands. */
  grid?: 'none' | 'strong'
  /** Step between ticks on a fixed value axis; absent shows its ends and middle. */
  valueStep?: number
  /** Fixed start and end of the value axis, so charts side by side share a scale. */
  valueRange?: [number, number]
  /** Label of a bracket joining the first two categories, such as a test's p-value. */
  compareLabel?: string
  /** A faint band behind each category whose label is named in the figure's colours,
   *  across the whole pane, so bands line up with the same rows beside them. */
  categoryBands?: boolean
  /** Leave out a grid chart's row headers, for labelling the rows beside it. */
  hideRowHeaders?: boolean
  /** Space between the columns of a grid chart, as a fraction of the pane. */
  gridGap?: number
  /** Space between its rows; room for the upper row's tick labels. */
  gridRowGap?: number
}

export type LegendPlacement = 'right' | 'top' | 'bottom' | 'each' | 'none'

export interface FigureGroup {
  id: string
  /** Letter drawn at the group's top left; null uses the next letter in order. */
  label: string | null
  rows: number
  cols: number
  /** Height relative to the other groups. */
  weight: number
  /** Width of each column relative to the others; absent gives equal widths. */
  colWeights?: number[]
  /** Height of each row relative to the others; absent gives equal heights. */
  rowWeights?: number[]
  /** CSS colour behind the group, or null for none. */
  background: string | null
  /** right: Plotly's legend on the top-right pane. top and bottom: one
   *  horizontal legend across the group. 'shared' is the old name for right. */
  legend: LegendPlacement | 'shared'
  /** Title drawn above the group, and above a top legend. */
  title?: string | null
  /** Sets the one-row legend in italics, as for genus and species names. */
  legendItalic?: boolean
  /** Sets the one-row legend in bold. */
  legendBold?: boolean
  panes: FigurePane[]
}

export type PageSize = 'a4' | 'a4-landscape' | 'letter' | 'letter-landscape' | 'single' | 'onehalf' | 'double' | 'custom'

export interface FigureDoc {
  id: string
  title: string
  page: {
    size: PageSize
    /** Used by the journal widths and custom sizes. */
    widthMm: number
    heightMm: number
    marginMm: number
  }
  /** titleBold sets chart and group titles in bold. */
  font: { family: 'Arial' | 'Times New Roman'; sizePt: number; titlePt?: number; tickPt?: number; titleBold?: boolean }
  /** Colour for every trace with this legend name, across all panes. */
  colours?: Record<string, string>
  /** Text shown in place of a run, group or category label on the axes and grid headers. */
  labels?: Record<string, string>
  /** none draws only the letters a group sets itself. */
  letters: { case: 'lower' | 'upper' | 'none'; sizePt: number }
  groupGapMm: number
  /** Resolution of the PNG saved beside the layout and of raster exports. */
  rasterDpi?: number
  groups: FigureGroup[]
}

export const PAGE_SIZES: Record<PageSize, { label: string; w: number; h: number; fixed: boolean }> = {
  'a4':               { label: 'A4 portrait',            w: 210,   h: 297,   fixed: true },
  'a4-landscape':     { label: 'A4 landscape',           w: 297,   h: 210,   fixed: true },
  'letter':           { label: 'Letter portrait',        w: 215.9, h: 279.4, fixed: true },
  'letter-landscape': { label: 'Letter landscape',       w: 279.4, h: 215.9, fixed: true },
  'single':           { label: 'Single column (85 mm)',  w: 85,    h: 100,   fixed: false },
  'onehalf':          { label: '1.5 column (120 mm)',    w: 120,   h: 140,   fixed: false },
  'double':           { label: 'Double column (180 mm)', w: 180,   h: 200,   fixed: false },
  'custom':           { label: 'Custom',                 w: 180,   h: 200,   fixed: false },
}

/** Page width and height in mm. Journal widths keep their width and take the chosen height. */
export function pageDims(doc: FigureDoc): { w: number; h: number } {
  const p = PAGE_SIZES[doc.page.size]
  if (p.fixed) return { w: p.w, h: p.h }
  return { w: doc.page.size === 'custom' ? doc.page.widthMm : p.w, h: doc.page.heightMm }
}

export const emptyPane = (): FigurePane => ({ chart: null, title: null, xTitle: null })

const newId = () => Math.random().toString(36).slice(2, 10)

export function newGroup(rows = 1, cols = 2): FigureGroup {
  return {
    id: newId(), label: null, rows, cols, weight: rows, background: null, legend: 'right',
    panes: Array.from({ length: rows * cols }, emptyPane),
  }
}

/** Resize a group's grid, keeping panes that still fit at their row and column. */
export function resizeGroup(g: FigureGroup, rows: number, cols: number): FigureGroup {
  const panes: FigurePane[] = []
  for (let r = 0; r < rows; r++)
    for (let c = 0; c < cols; c++)
      panes.push(r < g.rows && c < g.cols ? g.panes[r * g.cols + c] : emptyPane())
  const colWeights = g.colWeights && Array.from({ length: cols }, (_, c) => g.colWeights![c] ?? 1)
  const rowWeights = g.rowWeights && Array.from({ length: rows }, (_, r) => g.rowWeights![r] ?? 1)
  return { ...g, rows, cols, weight: g.weight === g.rows ? rows : g.weight, panes, colWeights, rowWeights }
}

/** Each column's share of the group width, summing to 1. */
export function columnShares(g: FigureGroup): number[] {
  const w = Array.from({ length: g.cols }, (_, c) => Math.max(0.1, g.colWeights?.[c] ?? 1))
  const total = w.reduce((a, b) => a + b, 0)
  return w.map(x => x / total)
}

export function newFigure(title: string): Omit<FigureDoc, 'id'> {
  return {
    title,
    page: { size: 'a4', widthMm: 180, heightMm: 200, marginMm: 10 },
    font: { family: 'Arial', sizePt: 8 },
    letters: { case: 'lower', sizePt: 16 },
    groupGapMm: 4,
    groups: [newGroup(1, 2)],
  }
}

/** The letter each group shows, counting only groups without their own label.
 *  Empty when letters are off and the group sets none. */
export function groupLetters(doc: FigureDoc): string[] {
  let n = 0
  return doc.groups.map(g => {
    if (g.label != null) return g.label
    if (doc.letters.case === 'none') return ''
    const l = String.fromCharCode(97 + n++)
    return doc.letters.case === 'upper' ? l.toUpperCase() : l
  })
}

export const MM_PER_IN = 25.4
/** CSS pixels per mm; charts are laid out at 96 px per inch. */
export const PX_PER_MM = 96 / MM_PER_IN
export const PX_PER_PT = 96 / 72

export interface Box { x: number; y: number; w: number; h: number }

export const legendPlacement = (g: FigureGroup): LegendPlacement => g.legend === 'shared' ? 'right' : g.legend

export const titlePt = (doc: FigureDoc) => doc.font.titlePt ?? doc.font.sizePt * 1.15
export const tickPt = (doc: FigureDoc) => doc.font.tickPt ?? doc.font.sizePt * 0.9

const PT_TO_MM = MM_PER_IN / 72

export interface TextRun { text: string; bold: boolean; italic: boolean }
export interface TextLineLayout {
  /** Baseline, in px from the pane's top. */
  y: number
  x: number
  anchor: 'start' | 'middle' | 'end'
  sizePx: number
  fill: string
  runs: TextRun[]
}

/** Split *italic*, **bold** and ***bold italic*** markup into runs. Unpaired markers stay as text. */
export function parseRuns(text: string): TextRun[] {
  const runs: TextRun[] = []
  const re = /\*\*\*(.+?)\*\*\*|\*\*(.+?)\*\*|\*(.+?)\*/g
  let at = 0
  for (let m = re.exec(text); m; m = re.exec(text)) {
    if (m.index > at) runs.push({ text: text.slice(at, m.index), bold: false, italic: false })
    runs.push(m[1] != null ? { text: m[1], bold: true, italic: true }
            : m[2] != null ? { text: m[2], bold: true, italic: false } : { text: m[3], bold: false, italic: true })
    at = re.lastIndex
  }
  if (at < text.length || !runs.length) runs.push({ text: text.slice(at), bold: false, italic: false })
  return runs
}

/** A picture pane's split: the picture's height, then its caption's, in the units of `h`.
 *  The caption takes its lines' natural height, at most half the pane. */
export function captionSplit(pane: FigurePane, doc: FigureDoc, h: number): { imageH: number; captionH: number } {
  if (!pane.image || !pane.text?.lines.length) return { imageH: h, captionH: 0 }
  const natural = pane.text.lines.reduce((a, l) => a + (l.sizePt ?? doc.font.sizePt) * PX_PER_PT * 1.2, 0)
  const captionH = Math.min(h / 2, natural)
  return { imageH: h - captionH, captionH }
}

/** Where each line of a text pane goes in a pane w × h px. Lines are 1.2 em apart.
 *  With `measure` (a run's width in px at a size), when a line is wider than 94% of the pane,
 *  every line shrinks by the same factor, so their sizes keep their proportions. */
export function layoutText(t: PaneText, doc: FigureDoc, w: number, h: number,
                           measure?: (run: TextRun, sizePx: number) => number): TextLineLayout[] {
  const runs = t.lines.map(l => parseRuns(l.text))
  const natural = t.lines.map(l => (l.sizePt ?? doc.font.sizePt) * PX_PER_PT)
  // Shrunk text keeps a small margin, so neighbouring panes do not touch.
  const fit = Math.min(1, ...natural.map((size, i) => {
    const width = measure ? runs[i].reduce((a, r) => a + measure(r, size), 0) : 0
    return width > w * 0.94 ? w * 0.94 / width : 1
  }))
  const sizes = natural.map(s => s * fit)
  const total = sizes.reduce((a, s) => a + s * 1.2, 0)
  let top = t.valign === 'top' ? 0 : t.valign === 'bottom' ? h - total : (h - total) / 2
  const anchor = t.align === 'left' ? 'start' : t.align === 'right' ? 'end' : 'middle'
  const x = t.align === 'left' ? 0 : t.align === 'right' ? w : w / 2
  return t.lines.map((l, i) => {
    const size = sizes[i]
    // The baseline sits where the line's 1.2 em box puts a cap height of about 0.7 em in the middle.
    const y = top + size * 0.95
    top += size * 1.2
    return { y, x, anchor, sizePx: size, fill: l.colour ?? '#222222', runs: runs[i] }
  })
}

export interface GroupGeometry {
  box: Box
  letter: string
  /** The group title's strip, when it has one. */
  title: Box | null
  /** The strip a top or bottom legend is drawn in. */
  legend: Box | null
  panes: Box[]
}

export interface FigureGeometry {
  page: { w: number; h: number }
  groups: GroupGeometry[]
}

/** Every group's and pane's box on the page, in mm. */
export function figureGeometry(doc: FigureDoc): FigureGeometry {
  const page = pageDims(doc)
  const m = doc.page.marginMm
  const innerW = page.w - 2 * m
  const innerH = page.h - 2 * m
  const gaps = Math.max(0, doc.groups.length - 1) * doc.groupGapMm
  const total = doc.groups.reduce((s, g) => s + Math.max(0.1, g.weight), 0) || 1
  const letters = groupLetters(doc)
  let y = m
  return {
    page,
    groups: doc.groups.map((g, i) => {
      const h = (innerH - gaps) * Math.max(0.1, g.weight) / total
      const box = { x: m, y, w: innerW, h }
      y += h + doc.groupGapMm
      let top = box.y, bottom = box.y + box.h
      let title: Box | null = null, legend: Box | null = null
      if (g.title) {
        const th = titlePt(doc) * 1.8 * PT_TO_MM
        title = { x: box.x, y: top, w: box.w, h: th }
        top += th
      }
      const place = legendPlacement(g)
      if (place === 'top' || place === 'bottom') {
        const lh = doc.font.sizePt * 2 * PT_TO_MM
        legend = { x: box.x, y: place === 'top' ? top : bottom - lh, w: box.w, h: lh }
        if (place === 'top') top += lh
        else bottom -= lh
      }
      const shares = columnShares(g)
      const lefts = shares.map((_, c) => shares.slice(0, c).reduce((a, b) => a + b, 0))
      const rw = Array.from({ length: g.rows }, (_, r) => Math.max(0.1, g.rowWeights?.[r] ?? 1))
      const rTotal = rw.reduce((a, b) => a + b, 0)
      const rh = rw.map(v => (bottom - top) * v / rTotal)
      const tops = rh.map((_, r) => top + rh.slice(0, r).reduce((a, b) => a + b, 0))
      const panes = g.panes.map((_, k) => ({
        x: box.x + lefts[k % g.cols] * box.w, y: tops[Math.floor(k / g.cols)],
        w: shares[k % g.cols] * box.w, h: rh[Math.floor(k / g.cols)],
      }))
      return { box, letter: letters[i], title, legend, panes }
    }),
  }
}
