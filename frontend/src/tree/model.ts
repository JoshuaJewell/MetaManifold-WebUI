// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).

// Tree topology, stable identifiers and the view-state operation layer.
//
// The source tree is parsed once and never changed. Every view (rerooted,
// collapsed, renamed) is rebuilt from it by replaying TreeOps. Branches are
// identified by their bipartition of the tip set, which does not depend on
// where the tree is rooted, so an op recorded under one rooting still finds
// its branch after a reroot or after the file is regenerated with the same
// reference taxa.

/** A node of the tree as written in the file. */
export interface SrcNode {
  name:     string
  length:   number | null
  /** jplace `{n}` edge number of the branch above this node. */
  edgeNum:  number | null
  support:  string | null
  children: number[]
  parent:   number
}

export interface Placement {
  edgeNum:  number
  lwr:      number
  distal:   number
  pendant:  number
  logLik:   number | null
}

export interface PQuery {
  names:  string[]
  /** Sequences this query stands for: the jplace `nm` multiplicities, else one per name. */
  mult:   number
  /** Candidates, best first. */
  places: Placement[]
}

export interface SourceTree {
  nodes:      SrcNode[]
  root:       number
  hasLengths: boolean
  pqueries:   PQuery[]
  meta:       Record<string, unknown>
}

export type TreeOp =
  | { op: 'reroot'; split: string | null; at: number }
  | { op: 'collapse'; clade: string }
  | { op: 'expand'; clade: string }
  | { op: 'expand_all' }
  /** Toggles the reversed order of a node's children. */
  | { op: 'rotate'; clade: string }
  /** `side` is the clade (tip-set hash) a branch name was given on; the root halves share one branch. */
  | { op: 'rename'; key: string; name: string | null; side?: string }
  | { op: 'style'; clade: string; style: NodeStyle | null }
  /** Edits copied from another tree's view, undone as one step. */
  | { op: 'batch'; from: string; ops: TreeOp[] }
  /** Support values imported from another tree, by split id. Replaces any earlier import. */
  | { op: 'support'; from: string | null; values: [string, string][] }

export interface Rgba { hex: string; alpha: number }

export interface NodeStyle {
  bold?:    boolean
  italic?:  boolean
  label?:   Rgba
  branch?:  Rgba
  /** Fill of the triangle drawn when the clade is collapsed. */
  fill?:    Rgba
  /** Descendants take this style unless they set their own. */
  inherit:  boolean
}

export type ResolvedStyle = Omit<NodeStyle, 'inherit'>

export interface ViewState {
  root:      { split: string; at: number } | null
  collapsed: Set<string>
  rotated:   Set<string>
  names:     Map<string, string>
  nameSides: Map<string, string>
  styles:    Map<string, NodeStyle>
  support:   Map<string, string>
  supportFrom: string | null
}

const SUPPORT_RE = /^\d+(\.\d+)?(\/\d+(\.\d+)?)*$/

export function parseNewick(text: string): SourceTree {
  const nodes: SrcNode[] = []
  const s = text.trim()
  let i = 0
  const mk = (parent: number): number => {
    nodes.push({ name: '', length: null, edgeNum: null, support: null, children: [], parent })
    return nodes.length - 1
  }
  const fail = (msg: string): never => { throw new Error(`Newick, character ${i}: ${msg}`) }
  const skipSpace = () => { while (i < s.length && /\s/.test(s[i])) i++ }

  // Bracket comments: a bare number is RAxML-style support, anything else is ignored.
  const readComment = (): string => {
    const end = s.indexOf(']', i)
    if (end < 0) fail('unclosed [')
    const body = s.slice(i + 1, end)
    i = end + 1
    return body
  }

  const readLabel = (): string => {
    skipSpace()
    if (s[i] === "'" || s[i] === '"') {
      const q = s[i++]
      let out = ''
      while (i < s.length) {
        if (s[i] === q) {
          if (s[i + 1] === q) { out += q; i += 2; continue }
          i++
          return out
        }
        out += s[i++]
      }
      fail('unclosed quote')
    }
    const start = i
    while (i < s.length && !'():,;[{'.includes(s[i])) i++
    return s.slice(start, i).trim()
  }

  const readTail = (n: number) => {
    const node = nodes[n]
    const label = readLabel()
    for (;;) {
      skipSpace()
      if (s[i] === ':') {
        i++
        skipSpace()
        const start = i
        while (i < s.length && /[-+0-9.eE]/.test(s[i])) i++
        const v = Number(s.slice(start, i))
        if (s.slice(start, i) === '' || !Number.isFinite(v)) fail('bad branch length')
        node.length = v
      } else if (s[i] === '{') {
        const end = s.indexOf('}', i)
        if (end < 0) fail('unclosed {')
        node.edgeNum = Number(s.slice(i + 1, end))
        i = end + 1
      } else if (s[i] === '[') {
        const c = readComment().trim()
        const v = SUPPORT_RE.test(c) ? c
          : c.startsWith('&&NHX') ? /:B=([^:\]]+)/.exec(c)?.[1] ?? null
          : c.startsWith('&') ? /(?:^|[&,])support=([^,\]]+)/.exec(c)?.[1] ?? null
          : null
        if (v !== null && node.support === null) node.support = v
      } else break
    }
    if (node.children.length > 0 && label !== '') {
      if (SUPPORT_RE.test(label)) node.support = label
      else node.name = label
    } else node.name = label
  }

  const root = mk(-1)
  const stack: number[] = []
  let cur = root
  skipSpace()
  // Leading comments such as FigTree's [&R] come before the tree.
  while (s[i] === '[') { readComment(); skipSpace() }
  if (s[i] !== '(') {
    readTail(root)
  } else {
    for (;;) {
      skipSpace()
      const c = s[i]
      if (c === '(') {
        i++
        stack.push(cur)
        const child = mk(cur)
        nodes[cur].children.push(child)
        cur = child
      } else if (c === ',') {
        i++
        const parent = stack[stack.length - 1]
        if (parent === undefined) fail('comma outside parentheses')
        const child = mk(parent)
        nodes[parent].children.push(child)
        cur = child
      } else if (c === ')') {
        i++
        const parent = stack.pop()
        if (parent === undefined) fail('unbalanced )')
        cur = parent!
        readTail(cur)
      } else if (c === ';' || i >= s.length) {
        break
      } else if (c === '[') {
        readComment()
      } else {
        const before = i
        readTail(cur)
        if (i === before) fail(`unexpected '${c}'`)
      }
    }
    if (stack.length) fail('unbalanced (')
  }
  const hasLengths = nodes.some((n, k) => k !== root && n.length !== null && n.length > 0)
  return { nodes, root, hasLengths, pqueries: [], meta: {} }
}

export function parseJplace(text: string): SourceTree {
  const doc = JSON.parse(text)
  const tree = parseNewick(String(doc.tree))
  const fields: string[] = doc.fields ?? []
  const col = (f: string) => fields.indexOf(f)
  const iEdge = col('edge_num'), iLwr = col('like_weight_ratio'), iDist = col('distal_length')
  const iPend = col('pendant_length'), iLik = col('likelihood')
  if (iEdge < 0) throw new Error('jplace fields lack edge_num')
  tree.pqueries = (doc.placements ?? []).map((pq: Record<string, unknown>) => {
    const rows = (pq.p ?? []) as number[][]
    const nm = Array.isArray(pq.nm) ? (pq.nm as [string, number][]) : null
    const names: string[] = Array.isArray(pq.n) ? (pq.n as string[]) : nm ? nm.map(x => x[0]) : []
    const mult = nm ? nm.reduce((a, x) => a + (Number(x[1]) || 0), 0) : Math.max(1, names.length)
    const places = rows.map(r => ({
      edgeNum: r[iEdge],
      lwr:     iLwr >= 0 ? r[iLwr] : 1,
      distal:  iDist >= 0 ? r[iDist] : 0,
      pendant: iPend >= 0 ? r[iPend] : 0,
      logLik:  iLik >= 0 ? r[iLik] : null,
    })).sort((a, b) => b.lwr - a.lwr)
    return { names, mult, places }
  })
  tree.meta = { ...(doc.metadata ?? {}), version: doc.version }
  return tree
}

export function parseTreeFile(format: 'jplace' | 'newick', text: string): SourceTree {
  return format === 'jplace' ? parseJplace(text) : parseNewick(text)
}

// Two independent 32-bit FNV-1a variants make a 64-bit tip hash. A tip set is
// hashed as the XOR of its tips, so a branch's two sides are h and total ^ h.
type H64 = [number, number]

function hashName(name: string): H64 {
  let a = 0x811c9dc5, b = 0x01000193 ^ 0x5bd1e995
  for (let k = 0; k < name.length; k++) {
    const c = name.charCodeAt(k)
    a = Math.imul(a ^ c, 0x01000193)
    b = Math.imul(b ^ c, 0x5bd1e995)
    b ^= b >>> 15
  }
  a ^= a >>> 13; a = Math.imul(a, 0x85ebca6b); a ^= a >>> 16
  b ^= b >>> 13; b = Math.imul(b, 0xc2b2ae35); b ^= b >>> 16
  return [a >>> 0, b >>> 0]
}

const hex = (h: H64) => h[0].toString(16).padStart(8, '0') + h[1].toString(16).padStart(8, '0')

/** Identifiers derived from the source tree. */
export interface TreeIndex {
  src:       SourceTree
  /** Unique key per tip: its name, suffixed `#2`, `#3` for repeats. */
  tipKey:    string[]
  tipCount:  number
  /** Hash of the tip set below each source node (tips on the child side of its branch). */
  side:      H64[]
  total:     H64
  /** Split id of the branch above each source node ('' for the root). */
  split:     string[]
  splitEdge: Map<string, number>
  /** Every tip set that sits on one side of some branch. */
  sides:     Set<string>
  edgeByNum: Map<number, number>
  /** True when the split id names the child side of the branch above this node. */
  childIsCanonical: boolean[]
  supportBySplit: Map<string, string>
}

export function indexTree(src: SourceTree): TreeIndex {
  const n = src.nodes.length
  const tipKey = new Array<string>(n).fill('')
  const seen = new Map<string, number>()
  let tipCount = 0
  const order: number[] = []
  const stack = [src.root]
  while (stack.length) {
    const v = stack.pop()!
    order.push(v)
    for (const c of src.nodes[v].children) stack.push(c)
  }
  for (const v of order) {
    if (src.nodes[v].children.length) continue
    tipCount++
    const name = src.nodes[v].name
    const k = (seen.get(name) ?? 0) + 1
    seen.set(name, k)
    tipKey[v] = k === 1 ? name : `${name}#${k}`
  }
  const refTip = order.filter(v => !src.nodes[v].children.length)
    .reduce((best, v) => (best < 0 || tipKey[v] < tipKey[best] ? v : best), -1)

  const side: H64[] = new Array(n)
  const hasRef = new Array<boolean>(n).fill(false)
  for (let k = order.length - 1; k >= 0; k--) {
    const v = order[k]
    const node = src.nodes[v]
    if (!node.children.length) {
      side[v] = hashName(tipKey[v])
      hasRef[v] = v === refTip
    } else {
      let a = 0, b = 0, r = false
      for (const c of node.children) { a ^= side[c][0]; b ^= side[c][1]; r ||= hasRef[c] }
      side[v] = [a >>> 0, b >>> 0]
      hasRef[v] = r
    }
  }
  const total = side[src.root]
  const split = new Array<string>(n).fill('')
  const childIsCanonical = new Array<boolean>(n).fill(false)
  const splitEdge = new Map<string, number>()
  const sides = new Set<string>()
  const edgeByNum = new Map<number, number>()
  const supportBySplit = new Map<string, string>()
  for (let v = 0; v < n; v++) {
    if (v === src.root) continue
    const other: H64 = [(total[0] ^ side[v][0]) >>> 0, (total[1] ^ side[v][1]) >>> 0]
    childIsCanonical[v] = !hasRef[v]
    split[v] = hex(hasRef[v] ? other : side[v])
    // A degree-two source root joins two branches with the same bipartition; the first wins.
    if (!splitEdge.has(split[v])) splitEdge.set(split[v], v)
    sides.add(hex(side[v]))
    sides.add(hex(other))
    const node = src.nodes[v]
    if (node.edgeNum !== null) edgeByNum.set(node.edgeNum, v)
    if (node.support !== null && node.children.length) supportBySplit.set(split[v], node.support)
  }
  return { src, tipKey, tipCount, side, total, split, splitEdge, sides, edgeByNum, childIsCanonical, supportBySplit }
}

/** Support values from another tree over the same tips, matched by bipartition. */
export function mapSupport(target: TreeIndex, supportTree: SourceTree): { support: Map<string, string>; matched: number; total: number } {
  const other = indexTree(supportTree)
  const support = new Map<string, string>()
  let total = 0, matched = 0
  for (const [s, v] of other.supportBySplit) {
    total++
    if (target.splitEdge.has(s)) { support.set(s, v); matched++ }
  }
  return { support, matched, total }
}

export function foldOps(ops: TreeOp[]): ViewState {
  const st: ViewState = {
    root: null, collapsed: new Set(), rotated: new Set(), names: new Map(), nameSides: new Map(), styles: new Map(), support: new Map(), supportFrom: null,
  }
  applyOps(st, ops)
  return st
}

function applyOps(st: ViewState, ops: TreeOp[]) {
  for (const o of ops) {
    switch (o.op) {
      case 'batch':      applyOps(st, o.ops); break
      case 'support':
        st.support = new Map(o.values)
        st.supportFrom = o.values.length ? o.from : null
        break
      case 'reroot':     st.root = o.split ? { split: o.split, at: o.at } : null; break
      case 'collapse':   st.collapsed.add(o.clade); break
      case 'expand':     st.collapsed.delete(o.clade); break
      case 'expand_all': st.collapsed.clear(); break
      case 'rotate':
        if (st.rotated.has(o.clade)) st.rotated.delete(o.clade)
        else st.rotated.add(o.clade)
        break
      case 'rename':
        if (o.name === null || o.name === '') { st.names.delete(o.key); st.nameSides.delete(o.key) }
        else {
          st.names.set(o.key, o.name)
          if (o.side) st.nameSides.set(o.key, o.side)
          else st.nameSides.delete(o.key)
        }
        break
      case 'style':
        if (o.style) st.styles.set(o.clade, o.style)
        else st.styles.delete(o.clade)
        break
    }
  }
}

/** The shortest op list that reproduces a view state. */
export function stateOps(st: ViewState): TreeOp[] {
  const ops: TreeOp[] = []
  if (st.root) ops.push({ op: 'reroot', split: st.root.split, at: st.root.at })
  for (const clade of st.collapsed) ops.push({ op: 'collapse', clade })
  for (const clade of st.rotated) ops.push({ op: 'rotate', clade })
  for (const [key, name] of st.names) {
    const side = st.nameSides.get(key)
    ops.push(side ? { op: 'rename', key, name, side } : { op: 'rename', key, name })
  }
  for (const [clade, style] of st.styles) ops.push({ op: 'style', clade, style })
  if (st.support.size) ops.push({ op: 'support', from: st.supportFrom, values: [...st.support] })
  return ops
}

/** Split another tree's view into the edits that apply to this tree and those that do not. */
export function matchOps(ix: TreeIndex, ops: TreeOp[]): { matched: TreeOp[]; unmatched: TreeOp[] } {
  const matched: TreeOp[] = [], unmatched: TreeOp[] = []
  for (const o of stateOps(foldOps(ops))) {
    if (o.op === 'support') {
      const values = o.values.filter(([s]) => ix.splitEdge.has(s))
      if (values.length) matched.push({ ...o, values })
      if (values.length < o.values.length) unmatched.push({ ...o, values: o.values.filter(([s]) => !ix.splitEdge.has(s)) })
      continue
    }
    (unmatchedOps(ix, [o]) ? unmatched : matched).push(o)
  }
  return { matched, unmatched }
}

/** Ops whose branch or tip does not exist in this tree. */
export function unmatchedOps(ix: TreeIndex, ops: TreeOp[]): number {
  const tips = new Set(ix.tipKey.filter(k => k !== ''))
  const flat = (list: TreeOp[]): TreeOp[] => list.flatMap(o => o.op === 'batch' ? flat(o.ops) : [o])
  return flat(ops).filter(o => {
    if (o.op === 'reroot') return o.split !== null && !ix.splitEdge.has(o.split)
    if (o.op === 'collapse' || o.op === 'expand' || o.op === 'style' || o.op === 'rotate')
      return !ix.sides.has(o.clade) && o.clade !== hex(ix.total)
    if (o.op === 'rename') {
      if (o.key.startsWith('tip:')) return !tips.has(o.key.slice(4))
      if (o.key.startsWith('split:')) return !ix.splitEdge.has(o.key.slice(6))
    }
    if (o.op === 'support') return o.values.some(([s]) => !ix.splitEdge.has(s))
    return false
  }).length
}

/** A node of the displayed (rooted) tree. */
export interface DNode {
  id:        number
  /** Source node, or -1 for a root placed on a branch. */
  src:       number
  children:  DNode[]
  parent:    DNode | null
  len:       number
  /** Source branch this node hangs from (a source node index), or -1 at the root. */
  edge:      number
  /** Distance from the source child end of `edge`, at this node and at its parent. */
  dChild:    number
  dParent:   number
  clade:     string
  split:     string
  tips:      number
  tipKey:    string | null
  name:      string
  support:   string | null
  collapsed: boolean
  style:     ResolvedStyle
  /** Style set on this node itself. */
  ownStyle:  NodeStyle | null
  /** Own length of the source branch segment this node was built from. */
  seg:       number
  /** Length of single-child nodes spliced into this branch above the segment. */
  above:     number
  /** Set on a spliced single-child node: the node whose branch now includes it. */
  mergedInto: DNode | null
}

export interface DisplayTree {
  root:  DNode
  /** Clades whose children are drawn in reverse order. */
  rotated: Set<string>
  nodes: DNode[]
  /** Display nodes carrying each source branch (two when the root sits on it). */
  byEdge: Map<number, DNode[]>
}

export function buildDisplay(ix: TreeIndex, st: ViewState, support: Map<string, string>): DisplayTree {
  const { src } = ix
  const nodes: DNode[] = []
  const byEdge = new Map<number, DNode[]>()
  const len = (v: number) => src.nodes[v].length ?? (src.hasLengths ? 0 : 1)

  const make = (srcIdx: number, parent: DNode | null, edge: number, l: number, dChild: number, dParent: number): DNode => {
    const d: DNode = {
      id: nodes.length, src: srcIdx, children: [], parent, len: l, edge, dChild, dParent,
      clade: '', split: edge >= 0 ? ix.split[edge] : '', tips: 0,
      tipKey: srcIdx >= 0 && !src.nodes[srcIdx].children.length ? ix.tipKey[srcIdx] : null,
      name: '', support: null, collapsed: false, style: {}, ownStyle: null, seg: l, above: 0, mergedInto: null,
    }
    nodes.push(d)
    if (edge >= 0) {
      const list = byEdge.get(edge)
      if (list) list.push(d)
      else byEdge.set(edge, [d])
    }
    return d
  }

  // Walk the unrooted graph away from `from`, entering source node `v` over branch `edge`.
  const grow = (v: number, from: number, parent: DNode, edge: number, l: number, dChild: number, dParent: number) => {
    const d = make(v, parent, edge, l, dChild, dParent)
    parent.children.push(d)
    const node = src.nodes[v]
    for (const c of node.children) if (c !== from) grow(c, v, d, c, len(c), 0, len(c))
    if (node.parent >= 0 && node.parent !== from) grow(node.parent, v, d, v, len(v), len(v), 0)
  }

  let root: DNode
  const rootEdge = st.root ? ix.splitEdge.get(st.root.split) : undefined
  if (st.root && rootEdge !== undefined) {
    const c = rootEdge, p = src.nodes[c].parent
    const L = len(c)
    // `at` is measured from the canonical side, which may be either end.
    const fromChild = ix.childIsCanonical[c] ? st.root.at * L : (1 - st.root.at) * L
    root = make(-1, null, -1, 0, 0, 0)
    grow(c, p, root, c, fromChild, 0, fromChild)
    grow(p, c, root, c, L - fromChild, L, fromChild)
  } else {
    root = make(src.root, null, -1, 0, 0, 0)
    for (const c of src.nodes[src.root].children) grow(c, src.root, root, c, len(c), 0, len(c))
  }

  // A degree-two vertex (a rooted file's root, once rerooted away) would draw as
  // a node with one child; fold it into the branch below.
  const splice = (d: DNode) => {
    for (let k = 0; k < d.children.length; k++) {
      let ch = d.children[k]
      while (ch.children.length === 1) {
        const g = ch.children[0]
        g.above = ch.above + ch.seg
        g.len += ch.len
        g.parent = d
        ch.mergedInto = g
        d.children[k] = g
        ch = g
      }
      splice(ch)
    }
  }
  splice(root)

  const post: DNode[] = []
  const walk = [root]
  while (walk.length) {
    const d = walk.pop()!
    post.push(d)
    for (const c of d.children) walk.push(c)
  }
  for (let k = post.length - 1; k >= 0; k--) {
    const d = post[k]
    if (d.tipKey !== null) {
      d.clade = hex(ix.side[d.src])
      d.tips = 1
      d.name = st.names.get(`tip:${d.tipKey}`) ?? src.nodes[d.src].name
    } else {
      let a = 0, b = 0, t = 0
      for (const c of d.children) {
        const h = parseClade(c.clade)
        a ^= h[0]; b ^= h[1]; t += c.tips
      }
      d.clade = hex([a >>> 0, b >>> 0])
      d.tips = t
      d.name = d.split ? st.names.get(`split:${d.split}`) ?? '' : ''
      d.support = d.split ? support.get(d.split) ?? null : null
      d.collapsed = d.parent !== null && st.collapsed.has(d.clade)
    }
  }

  // A root placed on a branch splits it between two children with the same split id.
  // Its name and support belong on one of them: the side the name was given on,
  // else the smaller side.
  const carriers = new Map<string, DNode[]>()
  for (const d of post) {
    if (d.tipKey !== null || !d.split) continue
    const list = carriers.get(d.split)
    if (list) list.push(d)
    else carriers.set(d.split, [d])
  }
  for (const [split, list] of carriers) {
    if (list.length < 2) continue
    const side = st.nameSides.get(`split:${split}`)
    const owner = list.find(d => d.clade === side) ?? list.reduce((a, b) => (b.tips < a.tips ? b : a))
    for (const d of list) if (d !== owner) { d.name = ''; d.support = null }
  }
  // Styles resolve top-down: a node's own style overrides what it inherits.
  const clean = (x: NodeStyle): ResolvedStyle => {
    const out: ResolvedStyle = {}
    if (x.bold !== undefined) out.bold = x.bold
    if (x.italic !== undefined) out.italic = x.italic
    if (x.label) out.label = x.label
    if (x.branch) out.branch = x.branch
    if (x.fill) out.fill = x.fill
    return out
  }
  const paint = (d: DNode, inherited: ResolvedStyle) => {
    const own = st.styles.get(d.clade) ?? null
    d.ownStyle = own
    d.style = own ? { ...inherited, ...clean(own) } : inherited
    const pass = own?.inherit ? d.style : inherited
    for (const c of d.children) paint(c, pass)
  }
  paint(root, {})

  return { root, nodes: post, byEdge, rotated: st.rotated }
}

/** Where a point `distal` from the source child end of `d`'s branch lies on the drawn branch. */
export function pointOnBranch(d: DNode, distal: number): { node: DNode; fromParent: number } {
  const offset = d.above + Math.min(d.seg, Math.abs(distal - d.dParent))
  let node = d
  while (node.mergedInto) node = node.mergedInto
  return { node, fromParent: offset }
}

function parseClade(h: string): H64 {
  return [parseInt(h.slice(0, 8), 16), parseInt(h.slice(8), 16)]
}

export function nodeKey(d: DNode): string | null {
  if (d.tipKey !== null) return `tip:${d.tipKey}`
  return d.split ? `split:${d.split}` : null
}

/** The reroot op that places the root halfway along the branch above `d`. */
export function rerootAbove(ix: TreeIndex, d: DNode): TreeOp | null {
  if (d.edge < 0) return null
  return { op: 'reroot', split: ix.split[d.edge], at: 0.5 }
}

/** Root at the midpoint of the longest tip-to-tip path. */
export function midpointOp(ix: TreeIndex): TreeOp | null {
  const { src } = ix
  const n = src.nodes.length
  const len = (v: number) => src.nodes[v].length ?? (src.hasLengths ? 0 : 1)
  const neighbours = (v: number): [number, number][] => {
    const out: [number, number][] = src.nodes[v].children.map(c => [c, c])
    if (src.nodes[v].parent >= 0) out.push([src.nodes[v].parent, v])
    return out
  }
  const farthest = (start: number) => {
    const dist = new Float64Array(n).fill(-1)
    const via = new Int32Array(n).fill(-1)
    dist[start] = 0
    const stack = [start]
    while (stack.length) {
      const v = stack.pop()!
      for (const [w, e] of neighbours(v)) {
        if (dist[w] >= 0) continue
        dist[w] = dist[v] + len(e)
        via[w] = v
        stack.push(w)
      }
    }
    let best = start
    for (let v = 0; v < n; v++) if (!src.nodes[v].children.length && dist[v] > dist[best]) best = v
    return { best, dist, via }
  }
  const firstTip = ix.tipKey.findIndex(k => k !== '')
  if (firstTip < 0) return null
  const a = farthest(firstTip).best
  const { best: b, dist, via } = farthest(a)
  const half = dist[b] / 2
  let v = b
  while (via[v] >= 0 && dist[via[v]] >= half) v = via[v]
  const u = via[v]
  if (u < 0) return null
  // The branch between u and v holds the midpoint; its source child is whichever is lower.
  const child = src.nodes[v].parent === u ? v : u
  const L = len(child)
  if (L <= 0) return { op: 'reroot', split: ix.split[child], at: 0.5 }
  const fromV = dist[v] - half
  const fromChild = child === v ? fromV : L - fromV
  const at = ix.childIsCanonical[child] ? fromChild / L : 1 - fromChild / L
  return { op: 'reroot', split: ix.split[child], at: Math.min(1, Math.max(0, at)) }
}

/** Tip keys below a display node, in display order. */
export function tipsBelow(d: DNode): DNode[] {
  const out: DNode[] = []
  const stack = [d]
  while (stack.length) {
    const x = stack.pop()!
    if (x.tipKey !== null) out.push(x)
    for (let k = x.children.length - 1; k >= 0; k--) stack.push(x.children[k])
  }
  return out
}

function quote(name: string): string {
  if (name === '') return ''
  return /^[^\s():,;[\]'"{}]+$/.test(name) ? name : `'${name.replace(/'/g, "''")}'`
}

const fmtLen = (x: number) => Number(x.toPrecision(10)).toString()

/**
 * The displayed tree as Newick, with renames applied. Collapsed clades are written in full.
 * Internal labels hold support or names; 'both' writes the name and puts support in an NHX comment.
 */
export function toNewick(t: DisplayTree, internal: 'support' | 'names' | 'both', lengths: boolean): string {
  const write = (d: DNode): string => {
    let s = ''
    if (d.children.length) s = '(' + d.children.map(write).join(',') + ')'
    const label = d.tipKey !== null ? d.name
      : internal === 'support' ? d.support ?? '' : d.name
    s += quote(label)
    if (lengths && d.parent) s += ':' + fmtLen(d.len)
    if (internal === 'both' && d.tipKey === null && d.support !== null) s += `[&&NHX:B=${d.support}]`
    return s
  }
  return write(t.root) + ';'
}
