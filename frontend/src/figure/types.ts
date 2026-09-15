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

export interface FigurePane {
  chart: ChartSpec | null
  /** Replaces the chart's own title; empty hides it, null keeps it. */
  title: string | null
  /** Text for the x-axis title; null keeps the chart's own. */
  xTitle: string | null
  style?: PaneStyle
}

/** Adjustments to the chart's own look; an absent key keeps the chart's. */
export interface PaneStyle {
  /** x tick label angle in degrees. */
  tickAngle?: number
  /** Label every x category instead of letting Plotly skip some. */
  allTicks?: boolean
  /** Opacity of line-only traces, such as the lines joining paired samples. */
  lineOpacity?: number
  /** Outline width of sample points, in px. */
  pointBorder?: number
  /** Opacity of sample points. */
  pointOpacity?: number
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
  /** CSS colour behind the group, or null for none. */
  background: string | null
  /** right: Plotly's legend on the top-right pane. top and bottom: one
   *  horizontal legend across the group. 'shared' is the old name for right. */
  legend: LegendPlacement | 'shared'
  /** Title drawn above the group, and above a top legend. */
  title?: string | null
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
  font: { family: 'Arial' | 'Times New Roman'; sizePt: number; titlePt?: number; tickPt?: number }
  /** Colour for every trace with this legend name, across all panes. */
  colours?: Record<string, string>
  letters: { case: 'lower' | 'upper'; sizePt: number }
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
  return { ...g, rows, cols, weight: g.weight === g.rows ? rows : g.weight, panes }
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

/** The letter each group shows, counting only groups without their own label. */
export function groupLetters(doc: FigureDoc): string[] {
  let n = 0
  return doc.groups.map(g => {
    if (g.label != null) return g.label
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
      const pw = box.w / g.cols, ph = (bottom - top) / g.rows
      const panes = g.panes.map((_, k) => ({
        x: box.x + (k % g.cols) * pw, y: top + Math.floor(k / g.cols) * ph, w: pw, h: ph,
      }))
      return { box, letter: letters[i], title, legend, panes }
    }),
  }
}
