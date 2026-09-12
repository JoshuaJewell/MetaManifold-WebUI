// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { symbolPath, type Dataset, type TreeLayout, type TreeSettings } from './layout'
import { fmt } from './ui'

export interface LegendRow { dataset: Dataset; label: string; values: number[]; sizes: number[]; height: number }

export function legendItems(layout: TreeLayout, s: TreeSettings): LegendRow[] {
  if (!s.legend) return []
  const rows: LegendRow[] = []
  for (const ds of ['placements', 'support'] as Dataset[]) {
    const sc = layout.scales[ds]
    if (!sc || !s.symbols[ds].show || sc.max <= 0) continue
    const nice = (v: number) => Number(v.toPrecision(2))
    const raw = ds === 'placements' && s.placements === 'best'
      ? [1, Math.round(sc.max / 2), sc.max]
      : [nice(sc.max / 4), nice(sc.max / 2), sc.max]
    const values = [...new Set(raw.filter(v => v > 0))]
    const sizes = values.map(sc.size)
    rows.push({ dataset: ds, label: sc.label, values, sizes, height: Math.max(...sizes, s.fontSize) + 10 })
  }
  return rows
}

export function Legend({ rows, settings }: { rows: LegendRow[]; settings: TreeSettings }) {
  const fs = settings.fontSize
  let y = 8
  return (
    <g>
      {rows.map(row => {
        const st = settings.symbols[row.dataset]
        const mid = y + row.height / 2
        y += row.height
        let x = 20 + row.label.length * fs * 0.6 + 12
        return (
          <g key={row.dataset}>
            <text x={20} y={mid} dominantBaseline="central" fill="currentColor">{row.label}</text>
            {row.values.map((v, k) => {
              const size = row.sizes[k]
              const cxs = x + size / 2
              x += size + 8 + String(fmt(v, 3)).length * fs * 0.6 + 10
              return (
                <g key={k}>
                  <path d={symbolPath(st.shape, cxs, mid, size)} fill={st.fill.hex} fillOpacity={st.fill.alpha}
                    stroke={st.stroke.hex} strokeOpacity={st.stroke.alpha} strokeWidth={1} />
                  <text x={cxs + size / 2 + 4} y={mid} dominantBaseline="central" fill="currentColor" fontSize={fs * 0.9}>{fmt(v, 3)}</text>
                </g>
              )
            })}
          </g>
        )
      })}
    </g>
  )
}
