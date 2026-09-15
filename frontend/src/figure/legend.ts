// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import type { PlotSpec } from './charts'

export interface LegendEntry {
  name: string
  fill: string
  stroke: string | null
  shape: 'square' | 'circle'
}

const str = (v: unknown): string | null => typeof v === 'string' ? v : null

/** The entries a chart would list in its own legend, in trace order. */
export function legendEntries(fig: PlotSpec): LegendEntry[] {
  const out: LegendEntry[] = []
  for (const t of fig.data) {
    const name = str(t.name)
    if (!name || t.showlegend === false || t.visible === false) continue
    const marker = (t.marker ?? {}) as Record<string, unknown>
    const line = (t.line ?? {}) as Record<string, unknown>
    if (t.type === 'box') {
      out.push({ name, fill: str(t.fillcolor) ?? str(marker.color) ?? '#999', stroke: str(line.color), shape: 'square' })
    } else if (t.type === 'scatter') {
      if (!String(t.mode ?? 'markers').includes('markers')) continue
      out.push({ name, fill: str(marker.color) ?? '#999', stroke: null, shape: 'circle' })
    } else {
      out.push({ name, fill: str(marker.color) ?? '#999', stroke: null, shape: 'square' })
    }
  }
  return out
}

/** Entries from several charts, each name once, in first-seen order. */
export function mergeEntries(lists: LegendEntry[][]): LegendEntry[] {
  const seen = new Map<string, LegendEntry>()
  for (const l of lists) for (const e of l) if (!seen.has(e.name)) seen.set(e.name, e)
  return [...seen.values()]
}

export interface LegendLayout {
  fontPx: number
  swatch: number
  items: { entry: LegendEntry; x: number; textX: number }[]
}

let ctx: CanvasRenderingContext2D | null = null
export function textWidth(text: string, fontPx: number, family: string): number {
  ctx ??= document.createElement('canvas').getContext('2d')
  if (!ctx) return text.length * fontPx * 0.55
  ctx.font = `${fontPx}px ${family}`
  return ctx.measureText(text).width
}

/** One centred row across `width` px. Spacing, then swatches, then text shrink until it fits. */
export function layoutLegend(entries: LegendEntry[], width: number, fontPx: number, family: string): LegendLayout {
  const labels = entries.map(e => e.name.replace(/_/g, ' '))
  for (let scale = 1; scale > 0.4; scale -= 0.05) {
    const f = fontPx * Math.min(1, scale + 0.3)
    const swatch = f * Math.max(0.6, scale)
    const gap = f * 0.35
    const spacing = f * 1.4 * scale
    const widths = labels.map(l => swatch + gap + textWidth(l, f, family))
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
