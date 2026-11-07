// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

import type { PlotSpec } from './charts'
import type { FigurePane } from './types'

export interface LegendEntry {
  name: string
  fill: string
  stroke: string | null
  shape: 'square' | 'circle'
  /** Set in italics, as genus and species names are. */
  italic?: boolean
  /** Too small to see in this chart: it holds its place in the order but is not listed for it. */
  faint?: boolean
}

const str = (v: unknown): string | null => typeof v === 'string' ? v : null

const barValues = (t: Record<string, unknown>) =>
  ((t.orientation === 'h' ? t.x : t.y) as unknown[] | undefined ?? []).map(Number).filter(Number.isFinite)

/** The entries a chart would list in its own legend, in trace order. A bar series that
 *  never reaches 0.5% of the tallest stack is left out: it cannot be seen to be keyed. */
export function legendEntries(fig: PlotSpec): LegendEntry[] {
  const out: LegendEntry[] = []
  const stacks: number[] = []
  for (const t of fig.data) if (t.type === 'bar') barValues(t).forEach((v, i) => { stacks[i] = (stacks[i] ?? 0) + Math.abs(v) })
  const floor = 0.005 * Math.max(0, ...stacks)
  for (const t of fig.data) {
    const name = str(t.name)
    if (!name || t.showlegend === false || t.visible === false) continue
    const faint = t.type === 'bar' && barValues(t).every(v => Math.abs(v) < floor)
    const marker = (t.marker ?? {}) as Record<string, unknown>
    const line = (t.line ?? {}) as Record<string, unknown>
    if (t.type === 'box') {
      out.push({ name, fill: str(t.fillcolor) ?? str(marker.color) ?? '#999', stroke: str(line.color), shape: 'square' })
    } else if (t.type === 'scatter') {
      if (!String(t.mode ?? 'markers').includes('markers')) continue
      out.push({ name, fill: str(marker.color) ?? '#999', stroke: null, shape: 'circle' })
    } else {
      out.push({ name, fill: str(marker.color) ?? '#999', stroke: null, shape: 'square', ...(faint ? { faint } : {}) })
    }
  }
  return out
}

/** A pane's legend entries: none when it is left out of the legend, and genus and
 *  species names in italics. */
export function paneLegendEntries(fig: PlotSpec, pane: FigurePane): LegendEntry[] {
  if (pane.style?.hideFromLegend) return []
  const c = pane.chart
  const taxa = c?.kind === 'composition' && c.tag === 'rank' && /^(genus|species)$/i.test(c.value)
  return legendEntries(fig).map(e => taxa && e.name !== 'Other' ? { ...e, italic: true } : e)
}

/** Entries from several charts, each name once, in first-seen order. A name too small to
 *  see in one chart keeps that chart's place in the order, and is listed if any chart shows it. */
export function mergeEntries(lists: LegendEntry[][]): LegendEntry[] {
  const seen = new Map<string, LegendEntry>()
  for (const l of lists) for (const e of l) {
    const had = seen.get(e.name)
    if (!had || (had.faint && !e.faint)) seen.set(e.name, e)
  }
  return [...seen.values()].filter(e => !e.faint)
}

export interface LegendLayout {
  fontPx: number
  swatch: number
  items: { entry: LegendEntry; x: number; textX: number }[]
}

let ctx: CanvasRenderingContext2D | null = null
export function textWidth(text: string, fontPx: number, family: string, italic = false, bold = false): number {
  ctx ??= document.createElement('canvas').getContext('2d')
  if (!ctx) return text.length * fontPx * (bold ? 0.6 : 0.55)
  ctx.font = `${italic ? 'italic ' : ''}${bold ? 'bold ' : ''}${fontPx}px ${family}`
  return ctx.measureText(text).width
}

/** One centred row across `width` px. Spacing, then swatches, then text shrink until it fits. */
export function layoutLegend(entries: LegendEntry[], width: number, fontPx: number, family: string, italic = false, bold = false): LegendLayout {
  const labels = entries.map(e => e.name.replace(/_/g, ' '))
  for (let scale = 1; scale > 0.4; scale -= 0.05) {
    const f = fontPx * Math.min(1, scale + 0.3)
    const swatch = f * Math.max(0.6, scale)
    const gap = f * 0.35
    const spacing = f * 1.4 * scale
    const widths = labels.map((l, i) => swatch + gap + textWidth(l, f, family, italic || entries[i].italic, bold))
    const total = widths.reduce((a, b) => a + b, 0) + spacing * Math.max(0, entries.length - 1)
    if (total <= width || scale <= 0.45) {
      let x = Math.max(0, (width - total) / 2)
      const items = entries.map((entry, i) => {
        const item = { entry: { ...entry, name: labels[i] }, x, textX: x + swatch + gap }
        x += widths[i] + spacing
        return item
      })
      return { fontPx: f, swatch, items }
    }
  }
  return { fontPx, swatch: fontPx, items: [] }
}
