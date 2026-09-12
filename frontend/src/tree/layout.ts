// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { pointOnBranch, type DisplayTree, type DNode, type PQuery, type Rgba, type TreeIndex } from './model'

export type SymbolShape = 'circle' | 'square' | 'triangle' | 'diamond'
export type Dataset = 'placements' | 'support'

export interface SymbolStyle {
  show:    boolean
  shape:   SymbolShape
  fill:    Rgba
  stroke:  Rgba
  /** Symbol diameter in px at the smallest and largest magnitude. */
  minSize: number
  maxSize: number
}

export type LayoutMode = 'rectangular' | 'slanted' | 'circular' | 'unrooted'

export interface TreeSettings {
  layout:            LayoutMode
  branchLengths:     boolean
  /** Drawn branch length: as in the file, capped, square root, or all equal. */
  branchTransform:   'none' | 'cap' | 'sqrt' | 'equal'
  /** Cap for long branches, in branch-length units; null uses the 99th percentile. */
  capAt:             number | null
  /** Circular layout: degrees of circle used, and where it starts clockwise from the top. */
  arc:               number
  arcStart:          number
  scaleBar:          boolean
  /** Scale bar length in branch-length units; null picks a round number. */
  scaleBarLength:    number | null
  ladderize:         'none' | 'up' | 'down'
  /** Clade size for ladderizing: every tip below, or drawn rows (a collapsed clade is one). */
  ladderizeBy:       'tips' | 'rows'
  fontSize:          number
  rowHeight:         number
  /** Tree width (rectangular) or radius (circular) in px, before zoom. */
  size:              number
  zoom:              number
  /** Support values written as text at nodes. */
  showSupport:       boolean
  supportMin:        number
  showInternalNames: boolean
  /** Count each query at its best placement, or spread it over candidates by LWR. */
  placements:        'best' | 'all'
  placementMinLwr:   number
  placementLabels:   boolean
  symbols:           Record<Dataset, SymbolStyle>
  legend:            boolean
  /** Collapsed clade shape: equal sides at the deepest tip, or sides at the shallowest and deepest tips. */
  wedge:             'isosceles' | 'scalene'
  collapsedCounts:   boolean
  /** Tip labels beside each tip, or in one column past the deepest tip. */
  tipLabels:         'tips' | 'aligned'
  tipJustify:        'left' | 'right'
  /** Newick file in the same study whose support values are mapped by bipartition. */
  supportTree:       string | null
  newickInternal:    'support' | 'names' | 'both'
}

export const DEFAULT_SETTINGS: TreeSettings = {
  layout: 'rectangular', branchLengths: true, branchTransform: 'none', capAt: null,
  arc: 360, arcStart: 0, scaleBar: true, scaleBarLength: null, ladderize: 'none', ladderizeBy: 'tips',
  fontSize: 11, rowHeight: 14, size: 700, zoom: 1,
  showSupport: true, supportMin: 0, showInternalNames: true,
  placements: 'best', placementMinLwr: 0.1, placementLabels: false, wedge: 'isosceles',
  collapsedCounts: true, tipLabels: 'tips', tipJustify: 'left',
  symbols: {
    placements: { show: true, shape: 'circle', fill: { hex: '#e8590c', alpha: 0.6 }, stroke: { hex: '#e8590c', alpha: 1 }, minSize: 5, maxSize: 18 },
    support:    { show: false, shape: 'circle', fill: { hex: '#495057', alpha: 0.8 }, stroke: { hex: '#495057', alpha: 0 }, minSize: 2, maxSize: 8 },
  },
  legend: true,
  supportTree: null, newickInternal: 'support',
}

export interface LNode {
  d:        DNode
  /** Position in px. */
  x:        number
  y:        number
  /** Parent position, and the elbow where the branch meets the parent's vertical or arc. */
  px:       number
  py:       number
  ex:       number
  ey:       number
  angle:    number
  radius:   number
  /** Collapsed wedge corners, when collapsed. */
  wedge?:   [number, number][]
  /** Drawn shorter than its length by the cap. */
  capped:   boolean
  kids:     LNode[]
}

export interface LPlacement {
  query:  number
  rank:   number
  lwr:    number
  x:      number
  y:      number
  /** Placements inside a collapsed clade are drawn at the clade. */
  hidden: boolean
  node:   DNode
}

/** One summary symbol: a dataset's magnitude on one drawn branch or node. */
export interface LSymbol {
  dataset: Dataset
  node:    DNode
  x:       number
  y:       number
  value:   number
  size:    number
  queries: { query: number; lwr: number; rank: number }[]
}

export interface DatasetScale {
  max:    number
  /** What the magnitude is, for the legend. */
  label:  string
  /** Diameter from magnitude. */
  size:   (v: number) => number
}

export interface TreeLayout {
  nodes:      LNode[]
  byId:       Map<number, LNode>
  placements: LPlacement[]
  symbols:    LSymbol[]
  scales:     Partial<Record<Dataset, DatasetScale>>
  width:      number
  height:     number
  mode:       LayoutMode
  circular:   boolean
  scaleBar:   { x: number; y: number; px: number; value: number } | null
  /** Transformed lengths make branch-length units non-linear. */
  capValue:   number | null
  /** Aligned label column: x (rectangular) or radius (circular), and its width. */
  alignAt:    number
  labelWidth: number
  cx:         number
  cy:         number
  /** px per unit of branch length */
  scale:      number
}

const charW = 0.58

export function layoutTree(t: DisplayTree, ix: TreeIndex, pqueries: PQuery[], s: TreeSettings): TreeLayout {
  const mode: LayoutMode = s.layout
  const circular = mode === 'circular'
  const unrooted = mode === 'unrooted'
  const rowCount = new Map<number, number>()
  const rowsOf = (d: DNode): number => {
    const hit = rowCount.get(d.id)
    if (hit !== undefined) return hit
    const n = d.tipKey !== null || d.collapsed ? 1 : d.children.reduce((a, c) => a + rowsOf(c), 0)
    rowCount.set(d.id, n)
    return n
  }
  const size = (d: DNode) => s.ladderizeBy === 'rows' ? rowsOf(d) : d.tips
  const order = (d: DNode) => {
    const kids = s.ladderize === 'none' ? d.children
      : [...d.children].sort((a, b) => s.ladderize === 'up' ? size(a) - size(b) : size(b) - size(a))
    return t.rotated.has(d.clade) ? [...kids].reverse() : kids
  }

  // Drawn length of each branch.
  const transform = s.branchLengths ? s.branchTransform : 'none'
  let capValue: number | null = null
  if (transform === 'cap') {
    if (s.capAt !== null && s.capAt > 0) capValue = s.capAt
    else {
      const lens = t.nodes.filter(d => d.parent).map(d => d.len).sort((a, b) => a - b)
      capValue = lens.length ? lens[Math.min(lens.length - 1, Math.floor(lens.length * 0.99))] : null
    }
  }
  const drawn = (d: DNode): number => {
    if (!s.branchLengths) return 1
    switch (transform) {
      case 'cap':   return capValue !== null ? Math.min(d.len, capValue) : d.len
      case 'sqrt':  return Math.sqrt(Math.max(0, d.len))
      case 'equal': return 1
      default:      return d.len
    }
  }

  // Visible structure, rows and depths in branch-length units.
  interface Raw { d: DNode; kids: Raw[]; depth: number; row: number; extent: number; near: number }
  let rows = 0
  let maxDepth = 0
  let maxLabel = 0
  const owner = new Map<number, DNode>()
  const height = new Map<number, number>()
  const hOf = (d: DNode): number => {
    if (d.tipKey !== null || d.collapsed) { height.set(d.id, 0); return 0 }
    let h = 0
    for (const c of d.children) h = Math.max(h, hOf(c) + 1)
    height.set(d.id, h)
    return h
  }
  const unit = drawn
  let cladoMax = 0
  if (!s.branchLengths) cladoMax = hOf(t.root)

  const markHidden = (d: DNode, top: DNode) => {
    for (const c of d.children) { owner.set(c.id, top); markHidden(c, top) }
  }
  const cladeExtent = (d: DNode): number => {
    let m = 0
    for (const c of d.children) m = Math.max(m, unit(c) + cladeExtent(c))
    return m
  }
  const cladeNearest = (d: DNode): number => {
    if (!d.children.length) return 0
    let m = Infinity
    for (const c of d.children) m = Math.min(m, unit(c) + cladeNearest(c))
    return m
  }

  const build = (d: DNode, depth: number): Raw => {
    const r: Raw = { d, kids: [], depth, row: 0, extent: 0, near: 0 }
    if (d.tipKey !== null || d.collapsed) {
      r.row = rows++
      if (d.collapsed) {
        markHidden(d, d)
        const far = cladeExtent(d)
        r.extent = s.branchLengths ? far : 1
        r.near = s.wedge === 'isosceles' ? r.extent
          : s.branchLengths ? cladeNearest(d) : far > 0 ? cladeNearest(d) / far : 1
      }
      const label = d.collapsed ? collapsedLabel(d, s.collapsedCounts) : d.name
      maxLabel = Math.max(maxLabel, label.length)
      maxDepth = Math.max(maxDepth, depth + r.extent)
      return r
    }
    for (const c of order(d)) {
      const cd = s.branchLengths ? depth + drawn(c) : cladoMax - (height.get(c.id) ?? 0)
      r.kids.push(build(c, cd))
    }
    r.row = (r.kids[0].row + r.kids[r.kids.length - 1].row) / 2
    maxDepth = Math.max(maxDepth, depth)
    return r
  }
  const rootDepth = s.branchLengths ? 0 : cladoMax - (height.get(t.root.id) ?? 0)
  const raw = build(t.root, rootDepth)
  if (maxDepth <= 0) maxDepth = 1

  const labelWidth = maxLabel * s.fontSize * charW
  const aligned = s.tipLabels === 'aligned' && !unrooted
  const labelPx = labelWidth + 16 + (aligned ? 12 : 0)
  const barSpace = s.scaleBar ? 34 : 0
  const isCapped = (d: DNode) => capValue !== null && d.parent !== null && d.len > capValue
  const nodes: LNode[] = []
  const byId = new Map<number, LNode>()
  let width: number, heightPx: number, cx = 0, cy = 0, scale: number, alignAt: number

  if (mode === 'rectangular' || mode === 'slanted') {
    const W = s.size * s.zoom
    const rowH = s.rowHeight * s.zoom
    scale = W / maxDepth
    const pad = 20
    const place = (r: Raw, parent: LNode | null): LNode => {
      const x = pad + r.depth * scale
      const y = pad + r.row * rowH
      const px = parent ? parent.x : x
      const py = parent ? parent.y : y
      const n: LNode = {
        d: r.d, x, y, px, py, ex: px, ey: mode === 'slanted' ? py : y,
        angle: 0, radius: 0, kids: [], capped: isCapped(r.d),
      }
      if (r.d.collapsed) {
        const x2 = x + r.extent * scale
        const h = Math.min(rowH * 0.45, 6 + Math.log2(r.d.tips) * 2)
        n.wedge = [[x, y], [x + r.near * scale, y - h], [x2, y + h]]
      }
      nodes.push(n)
      byId.set(r.d.id, n)
      n.kids = r.kids.map(k => place(k, n))
      return n
    }
    place(raw, null)
    alignAt = pad + W + 8
    width = pad + W + labelPx + pad
    heightPx = pad * 2 + Math.max(0, rows - 1) * rowH + barSpace
  } else if (unrooted) {
    // Equal-angle layout: each node gets a share of the circle in proportion to its drawn tips.
    const leaves = new Map<Raw, number>()
    const count = (r: Raw): number => {
      const n = r.kids.length ? r.kids.reduce((a, k) => a + count(k), 0) : 1
      leaves.set(r, n)
      return n
    }
    const total = count(raw)
    interface U { r: Raw; x: number; y: number; a: number; parent: U | null; kids: U[] }
    const us: U[] = []
    const lay = (r: Raw, parent: U | null, x: number, y: number, a: number, a0: number): U => {
      const u: U = { r, x, y, a, parent, kids: [] }
      us.push(u)
      let from = a0
      for (const k of r.kids) {
        const w = (leaves.get(k)! / total) * Math.PI * 2
        const mid = from + w / 2
        const len = Math.max(0, k.depth - r.depth)
        u.kids.push(lay(k, u, x + len * Math.cos(mid), y + len * Math.sin(mid), mid, from))
        from += w
      }
      return u
    }
    const top = lay(raw, null, 0, 0, 0, -Math.PI / 2)
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity
    for (const u of us) {
      const reach = u.r.d.collapsed ? u.r.extent : 0
      const fx = u.x + reach * Math.cos(u.a), fy = u.y + reach * Math.sin(u.a)
      minX = Math.min(minX, u.x, fx); maxX = Math.max(maxX, u.x, fx)
      minY = Math.min(minY, u.y, fy); maxY = Math.max(maxY, u.y, fy)
    }
    const span = Math.max(maxX - minX, maxY - minY, 1e-9)
    scale = s.size * s.zoom / span
    const margin = labelPx + 24
    const X = (x: number) => margin + (x - minX) * scale
    const Y = (y: number) => margin + (y - minY) * scale
    const place = (u: U, parent: LNode | null): LNode => {
      const x = X(u.x), y = Y(u.y)
      const px = parent ? parent.x : x, py = parent ? parent.y : y
      const n: LNode = { d: u.r.d, x, y, px, py, ex: px, ey: py, angle: u.a, radius: 0, kids: [], capped: isCapped(u.r.d) }
      if (u.r.d.collapsed) {
        const dx = Math.cos(u.a), dy = Math.sin(u.a)
        const h = 6 + Math.log2(u.r.d.tips) * 2
        const far = u.r.extent * scale, near = u.r.near * scale
        n.wedge = [[x, y], [x + near * dx - h * dy, y + near * dy + h * dx], [x + far * dx + h * dy, y + far * dy - h * dx]]
      }
      nodes.push(n)
      byId.set(u.r.d.id, n)
      n.kids = u.kids.map(k => place(k, n))
      return n
    }
    place(top, null)
    alignAt = 0
    width = margin * 2 + (maxX - minX) * scale
    heightPx = margin * 2 + (maxY - minY) * scale + barSpace
  } else {
    const R = s.size * s.zoom / 2
    scale = R / maxDepth
    const margin = labelPx + 24
    cx = cy = R + margin
    const full = s.arc >= 360
    const arcRad = Math.min(360, Math.max(10, s.arc)) * Math.PI / 180
    const span = full ? Math.PI * 2 * (rows > 1 ? (rows - 1) / rows : 1) : arcRad
    const start = -Math.PI / 2 + s.arcStart * Math.PI / 180
    const ang = (row: number) => start + (rows > 1 ? row / (rows - 1) : 0) * span
    const pt = (r: number, a: number): [number, number] => [cx + r * Math.cos(a), cy + r * Math.sin(a)]
    const place = (r: Raw, parent: LNode | null): LNode => {
      const a = ang(r.row)
      const rad = r.depth * scale
      const [x, y] = pt(rad, a)
      const pr = parent ? parent.radius : rad
      const [px, py] = parent ? [parent.x, parent.y] : [x, y]
      const [ex, ey] = pt(pr, a)
      const n: LNode = { d: r.d, x, y, px, py, ex, ey, angle: a, radius: rad, kids: [], capped: isCapped(r.d) }
      if (r.d.collapsed) {
        const r2 = rad + r.extent * scale
        const da = Math.min(span / Math.max(rows, 1) * 0.45, 0.08)
        n.wedge = [pt(rad, a), pt(rad + r.near * scale, a - da), pt(r2, a + da)]
      }
      nodes.push(n)
      byId.set(r.d.id, n)
      n.kids = r.kids.map(k => place(k, n))
      return n
    }
    place(raw, null)
    alignAt = R + 8
    width = cx * 2
    heightPx = cx * 2 + barSpace
  }

  // Placements: the distal length is measured from the source child end of the branch.
  const placements: LPlacement[] = []
  if (pqueries.length) {
    pqueries.forEach((q, qi) => {
      q.places.forEach((p, rank) => {
        if (s.placements === 'best' && rank > 0) return
        if (rank > 0 && p.lwr < s.placementMinLwr) return
        const e = ix.edgeByNum.get(p.edgeNum)
        if (e === undefined) return
        const cands = t.byEdge.get(e) ?? []
        const hit = cands.find(c => p.distal >= Math.min(c.dChild, c.dParent) - 1e-12
                                  && p.distal <= Math.max(c.dChild, c.dParent) + 1e-12) ?? cands[0]
        if (!hit) return
        const { node: dn, fromParent } = pointOnBranch(hit, p.distal)
        const top = owner.get(dn.id)
        const shown = top ? byId.get(top.id) : byId.get(dn.id)
        if (!shown) return
        let x: number, y: number
        if (top) {
          const w = shown.wedge
          ;[x, y] = w ? [(w[0][0] + w[1][0] + w[2][0]) / 3, (w[0][1] + w[1][1] + w[2][1]) / 3] : [shown.x, shown.y]
        } else {
          const f = dn.len > 0 ? Math.min(1, fromParent / dn.len) : 0.5
          x = shown.ex + (shown.x - shown.ex) * f
          y = shown.ey + (shown.y - shown.ey) * f
        }
        placements.push({ query: qi, rank, lwr: p.lwr, x, y, hidden: !!top, node: top ?? dn })
      })
    })
  }

  const symbols: LSymbol[] = []
  const scales: TreeLayout['scales'] = {}

  // Placement counts per drawn branch, weighted by multiplicity (and by LWR over all candidates).
  if (s.symbols.placements.show && placements.length) {
    const groups = new Map<number, LSymbol>()
    for (const p of placements) {
      const w = pqueries[p.query].mult * (s.placements === 'all' ? p.lwr : 1)
      let g = groups.get(p.node.id)
      if (!g) {
        g = { dataset: 'placements', node: p.node, x: 0, y: 0, value: 0, size: 0, queries: [] }
        groups.set(p.node.id, g)
      }
      g.x += p.x * w
      g.y += p.y * w
      g.value += w
      g.queries.push({ query: p.query, lwr: p.lwr, rank: p.rank })
    }
    let max = 0
    for (const g of groups.values()) max = Math.max(max, g.value)
    const st = s.symbols.placements
    const size = (v: number) => st.minSize + (st.maxSize - st.minSize) * Math.sqrt(max > 0 ? Math.min(1, v / max) : 0)
    scales.placements = { max, size, label: s.placements === 'all' ? 'Placements (LWR-weighted)' : 'Placements' }
    for (const g of groups.values()) {
      if (g.value <= 0) continue
      g.x /= g.value
      g.y /= g.value
      g.size = size(g.value)
      symbols.push(g)
    }
  }

  // Support at each drawn internal node; the scale runs to the largest value in the tree.
  if (s.symbols.support.show) {
    let max = 0
    for (const d of t.nodes) if (d.support !== null) max = Math.max(max, supportValue(d.support))
    if (max > 0) {
      const st = s.symbols.support
      const size = (v: number) => st.minSize + (st.maxSize - st.minSize) * Math.min(1, Math.max(0, v / max))
      scales.support = { max, size, label: 'Support' }
      for (const n of nodes) {
        const d = n.d
        if (d.support === null || !d.parent) continue
        const v = supportValue(d.support)
        if (!(v >= s.supportMin)) continue
        symbols.push({ dataset: 'support', node: d, x: n.x, y: n.y, value: v, size: size(v), queries: [] })
      }
    }
  }

  // The scale bar is only meaningful while drawn length is linear in branch length.
  let scaleBar: TreeLayout['scaleBar'] = null
  if (s.scaleBar && s.branchLengths && (transform === 'none' || transform === 'cap') && scale > 0) {
    const value = s.scaleBarLength && s.scaleBarLength > 0 ? s.scaleBarLength : niceLength(maxDepth / 5)
    scaleBar = { x: 20, y: heightPx - 16, px: value * scale, value }
  }

  return {
    nodes, byId, placements, symbols, scales, width, height: heightPx, mode, circular,
    scaleBar, capValue, alignAt, labelWidth, cx, cy, scale,
  }
}

/** A round length (1, 2 or 5 times a power of ten) near `x`. */
export function niceLength(x: number): number {
  if (!(x > 0)) return 1
  const p = Math.pow(10, Math.floor(Math.log10(x)))
  const m = x / p
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p
}

export const supportValue = (s: string) => Number(s.split('/')[0])

/** SVG path for a symbol of diameter `size` centred on (x, y). */
export function symbolPath(shape: SymbolShape, x: number, y: number, size: number): string {
  const r = size / 2
  switch (shape) {
    case 'square':   return `M${x - r},${y - r}h${size}v${size}h${-size}Z`
    case 'diamond':  return `M${x},${y - r}L${x + r},${y}L${x},${y + r}L${x - r},${y}Z`
    case 'triangle': {
      const h = r * 1.15
      return `M${x},${y - h}L${x + h * 0.95},${y + h * 0.6}L${x - h * 0.95},${y + h * 0.6}Z`
    }
    default:         return `M${x - r},${y}a${r},${r} 0 1 0 ${size},0a${r},${r} 0 1 0 ${-size},0`
  }
}

export function collapsedLabel(d: DNode, counts = true): string {
  if (!counts) return d.name
  return d.name ? `${d.name} (${d.tips})` : `${d.tips} tips`
}
