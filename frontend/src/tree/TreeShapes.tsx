// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useMemo } from 'react'
import type { DNode, Rgba } from './model'
import { collapsedLabel, supportValue, type LNode, type TreeLayout, type TreeSettings } from './layout'
import styles from './Tree.module.css'
import { fmt, SELECT_COLOUR } from './ui'

export function TreeShapes({ layout, settings, selSubtree, matchSet, onSelect, onToggle, onRename }: {
  layout: TreeLayout
  settings: TreeSettings
  selSubtree: Set<number>
  matchSet: Set<number>
  onSelect: (d: DNode) => void
  onToggle: (d: DNode) => void
  onRename: (d: DNode, target?: Element) => void
}) {
  const { circular, cx, cy, mode } = layout
  const fs = settings.fontSize
  const aligned = settings.tipLabels === 'aligned' && mode !== 'unrooted'
  // Rectangular branches are horizontal; every other mode draws straight from the elbow.
  const branchPath = (n: LNode) => mode === 'rectangular' ? `M${n.px},${n.y}H${n.x}` : `M${n.ex},${n.ey}L${n.x},${n.y}`
  const right = aligned && settings.tipJustify === 'right'

  // Branches are batched into one path per colour.
  const { strokes, selected, leaders, caps } = useMemo(() => {
    const byColour = new Map<string, { rgba: Rgba | null; d: string }>()
    let sel = '', lead = ''
    let caps = ''
    const arc = (n: LNode) => {
      const first = n.kids[0], last = n.kids[n.kids.length - 1]
      if (mode === 'slanted' || mode === 'unrooted') return ''
      if (!circular) return `M${n.x},${first.y}V${last.y}`
      const r = n.radius
      if (r <= 0) return ''
      const a1 = first.angle, a2 = last.angle
      const large = a2 - a1 > Math.PI ? 1 : 0
      return `M${cx + r * Math.cos(a1)},${cy + r * Math.sin(a1)}A${r},${r} 0 ${large} 1 ${cx + r * Math.cos(a2)},${cy + r * Math.sin(a2)}`
    }
    for (const n of layout.nodes) {
      let p = ''
      if (n.d.parent) p += branchPath(n)
      if (n.kids.length) p += arc(n)
      if (n.capped) {
        // Two slashes across the middle of a branch drawn shorter than it is.
        const sx = mode === 'rectangular' ? n.px : n.ex, sy = mode === 'rectangular' ? n.y : n.ey
        const L = Math.hypot(n.x - sx, n.y - sy) || 1
        const ux = (n.x - sx) / L, uy = (n.y - sy) / L
        const mx = (sx + n.x) / 2, my = (sy + n.y) / 2
        const vx = 4 * (-uy + 0.5 * ux), vy = 4 * (ux + 0.5 * uy)
        for (const o of [-2, 2]) {
          const qx = mx + o * ux, qy = my + o * uy
          caps += `M${qx - vx},${qy - vy}L${qx + vx},${qy + vy}`
        }
      }
      const c = n.d.style.branch ?? null
      const key = c ? `${c.hex}/${c.alpha}` : ''
      const entry = byColour.get(key)
      if (entry) entry.d += p
      else byColour.set(key, { rgba: c, d: p })
      if (selSubtree.has(n.d.id)) sel += p
      if (aligned && !n.kids.length) {
        const end = n.wedge ? farthest(n.wedge, circular, cx, cy) : null
        if (circular) {
          const r0 = end ? end.r : n.radius
          if (layout.alignAt - 4 > r0) lead += `M${cx + r0 * Math.cos(n.angle)},${cy + r0 * Math.sin(n.angle)}L${cx + (layout.alignAt - 4) * Math.cos(n.angle)},${cy + (layout.alignAt - 4) * Math.sin(n.angle)}`
        } else {
          const x0 = end ? end.x : n.x
          if (layout.alignAt - 4 > x0) lead += `M${x0},${n.y}H${layout.alignAt - 4}`
        }
      }
    }
    return { strokes: [...byColour.values()], selected: sel, leaders: lead, caps }
  }, [layout, circular, cx, cy, selSubtree, aligned, mode])

  // Label anchor: beside the node, or in the aligned column.
  const labelPos = (n: LNode, dx: number, reach?: { x: number; r: number; px?: number; py?: number }) => {
    if (mode === 'unrooted') {
      const deg = n.angle * 180 / Math.PI
      const flip = Math.cos(n.angle) < 0
      const bx = reach?.px ?? n.x, by = reach?.py ?? n.y
      const x = bx + dx * Math.cos(n.angle), y = by + dx * Math.sin(n.angle)
      return { x, y, anchor: flip ? 'end' as const : 'start' as const, transform: `rotate(${flip ? deg + 180 : deg} ${x} ${y})` }
    }
    if (!circular) {
      if (aligned) {
        return right
          ? { x: layout.alignAt + layout.labelWidth, y: n.y, anchor: 'end' as const, transform: undefined }
          : { x: layout.alignAt, y: n.y, anchor: 'start' as const, transform: undefined }
      }
      return { x: (reach ? reach.x : n.x) + dx, y: n.y, anchor: 'start' as const, transform: undefined }
    }
    const deg = n.angle * 180 / Math.PI
    const flip = Math.cos(n.angle) < 0
    const r = aligned ? (right ? layout.alignAt + layout.labelWidth : layout.alignAt) : (reach ? reach.r : n.radius) + dx
    const x = cx + r * Math.cos(n.angle), y = cy + r * Math.sin(n.angle)
    const atEnd = right ? !flip : flip
    return { x, y, anchor: atEnd ? 'end' as const : 'start' as const, transform: `rotate(${flip ? deg + 180 : deg} ${x} ${y})` }
  }

  const textStyle = (d: DNode, boldDefault = false) => ({
    fontStyle:   d.style.italic ? 'italic' : undefined,
    fontWeight:  (d.style.bold ?? boldDefault) ? 700 : undefined,
    fill:        d.style.label?.hex ?? 'currentColor',
    fillOpacity: d.style.label?.alpha,
  })

  return (
    <g>
      {strokes.map((b, k) => (
        <path key={k} d={b.d} fill="none" stroke={b.rgba?.hex ?? 'currentColor'} strokeOpacity={b.rgba?.alpha} strokeWidth={1} />
      ))}
      {caps && <path d={caps} fill="none" stroke="currentColor" strokeWidth={1.2} />}
      {layout.scaleBar && (() => {
        const b = layout.scaleBar
        return (
          <g>
            <path d={`M${b.x},${b.y - 4}V${b.y + 4}M${b.x},${b.y}H${b.x + b.px}M${b.x + b.px},${b.y - 4}V${b.y + 4}`}
              fill="none" stroke="currentColor" strokeWidth={1.2} />
            <text x={b.x + b.px / 2} y={b.y - 6} textAnchor="middle" fontSize={fs * 0.9} fill="currentColor">{fmt(b.value, 3)}</text>
          </g>
        )
      })()}
      {leaders && <path d={leaders} fill="none" stroke="currentColor" strokeOpacity={0.35} strokeWidth={0.6} strokeDasharray="1 2" />}
      {selected && <path data-ui="" d={selected} fill="none" stroke={SELECT_COLOUR} strokeWidth={2.5} />}
      {layout.nodes.map(n => {
        const d = n.d
        const hit = d.parent ? branchPath(n) : ''
        const isSel = selSubtree.has(d.id)
        const out: React.ReactNode[] = []
        if (hit) out.push(
          <path key="h" data-ui="" d={hit} stroke="transparent" strokeWidth={9} fill="none" className={styles.hit}
            onClick={e => { e.stopPropagation(); onSelect(d) }}
            onDoubleClick={e => { e.stopPropagation(); if (d.collapsed) onToggle(d); else onRename(d) }} />,
        )
        if (n.wedge) {
          const lp = labelPos(n, 4, farthest(n.wedge, circular, cx, cy))
          const label = collapsedLabel(d, settings.collapsedCounts)
          const fillC = isSel ? { hex: SELECT_COLOUR, alpha: 0.25 } : d.style.fill ?? (d.style.branch
            ? { hex: d.style.branch.hex, alpha: d.style.branch.alpha * 0.3 } : { hex: 'currentColor', alpha: 0.18 })
          out.push(
            <polygon key="w" points={n.wedge.map(p => p.join(',')).join(' ')}
              fill={fillC.hex} fillOpacity={fillC.alpha}
              stroke={d.style.branch?.hex ?? 'currentColor'} strokeOpacity={d.style.branch?.alpha} strokeWidth={1}
              className={styles.hit}
              onClick={e => { e.stopPropagation(); onSelect(d) }}
              onDoubleClick={e => { e.stopPropagation(); onToggle(d) }} />,
          )
          if (label) out.push(
            <text key="wl" x={lp.x} y={lp.y} textAnchor={lp.anchor} transform={lp.transform} dominantBaseline="central"
              {...textStyle(d)} className={styles.label}
              onClick={() => onSelect(d)}
              onDoubleClick={e => onRename(d, e.currentTarget)}>{label}</text>,
          )
        } else if (d.tipKey !== null) {
          const lp = labelPos(n, 4)
          const hitMatch = matchSet.has(d.id)
          out.push(
            <text key="t" x={lp.x} y={lp.y} textAnchor={lp.anchor} transform={lp.transform} dominantBaseline="central"
              {...textStyle(d)}
              {...(hitMatch ? { fill: '#e8590c', fillOpacity: 1, fontWeight: 700 } : {})}
              className={styles.label}
              onClick={() => onSelect(d)}
              onDoubleClick={e => onRename(d, e.currentTarget)}>{d.name}</text>,
          )
        } else {
          if (settings.showSupport && d.support !== null && d.parent && supportValue(d.support) >= settings.supportMin) {
            const x = circular ? n.x : n.x - 3
            const y = circular ? n.y : n.y - 3
            out.push(
              <text key="s" x={x} y={y} textAnchor={circular ? 'middle' : 'end'} fontSize={fs * 0.75}
                fill="currentColor" fillOpacity={0.7}>{d.support}</text>,
            )
          }
          if (settings.showInternalNames && d.name) {
            out.push(
              <text key="n" x={n.x + 3} y={circular ? n.y - 3 : n.y + fs * 0.9}
                fontSize={fs * 0.85} {...textStyle(d, true)} className={styles.label}
                onClick={() => onSelect(d)}
                onDoubleClick={e => onRename(d, e.currentTarget)}>{d.name}</text>,
            )
          }
        }
        return out.length ? <g key={d.id}>{out}</g> : null
      })}
    </g>
  )
}

/** The wedge corner farthest from the node: its x (rectangular), radius (circular) and point (unrooted). */
function farthest(w: [number, number][], circular: boolean, cx: number, cy: number) {
  const r = Math.max(Math.hypot(w[1][0] - cx, w[1][1] - cy), Math.hypot(w[2][0] - cx, w[2][1] - cy))
  const d1 = Math.hypot(w[1][0] - w[0][0], w[1][1] - w[0][1]), d2 = Math.hypot(w[2][0] - w[0][0], w[2][1] - w[0][1])
  const far = d1 > d2 ? w[1] : w[2]
  const ux = d1 || d2 ? (far[0] - w[0][0]) / Math.max(d1, d2) : 0, uy = d1 || d2 ? (far[1] - w[0][1]) / Math.max(d1, d2) : 0
  const reach = Math.max(d1, d2)
  return { x: Math.max(w[1][0], w[2][0]), r: circular ? r : 0, px: w[0][0] + ux * reach, py: w[0][1] + uy * reach }
}
