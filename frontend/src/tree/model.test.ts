// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { describe, expect, test } from 'bun:test'
import {
  buildDisplay, foldOps, indexTree, matchOps, midpointOp, parseJplace, parseNewick,
  pointOnBranch, stateOps, toNewick, type SourceTree, type TreeOp, type ViewState,
} from './model'

const splitSet = (src: SourceTree) => new Set(indexTree(src).splitEdge.keys())
const byName = (src: SourceTree, name: string) => src.nodes.findIndex(n => n.name === name)
const stateJson = (st: ViewState) => JSON.stringify({
  root: st.root,
  collapsed: [...st.collapsed].sort(),
  rotated: [...st.rotated].sort(),
  names: [...st.names].sort(),
  nameSides: [...st.nameSides].sort(),
  styles: [...st.styles].sort(),
  support: [...st.support].sort(),
})

describe('parseNewick', () => {
  test('quoted labels keep spaces and doubled quotes', () => {
    const t = parseNewick("('A b':1,'it''s':2,\"C d\":3);")
    const tips = t.nodes.filter(n => !n.children.length).map(n => n.name)
    expect(tips).toEqual(['A b', "it's", 'C d'])
  })

  test('a leading [&R] comment is skipped', () => {
    const t = parseNewick('[&R] ((A:1,B:1):1,C:2);')
    expect(t.nodes.filter(n => !n.children.length).map(n => n.name)).toEqual(['A', 'B', 'C'])
  })

  test('numeric internal labels are support, other labels are names', () => {
    const t = parseNewick('((A,B)95:1,(C,D)Clade:1,E);')
    const internal = t.nodes.filter((n, i) => n.children.length && i !== t.root)
    expect(internal.map(n => [n.name, n.support])).toEqual([['', '95'], ['Clade', null]])
  })

  test('NHX B= and bare bracket numbers are read as support', () => {
    const t = parseNewick('((A,B):1[&&NHX:B=88],(C,D):1[70],E);')
    const internal = t.nodes.filter((n, i) => n.children.length && i !== t.root)
    expect(internal.map(n => n.support)).toEqual(['88', '70'])
  })

  test('unbalanced input throws', () => {
    expect(() => parseNewick('((A,B);')).toThrow()
  })
})

describe('parseJplace', () => {
  test('edge numbers, fields and multiplicities', () => {
    const doc = {
      version: 3,
      fields: ['edge_num', 'like_weight_ratio', 'distal_length', 'pendant_length'],
      tree: '((A:1{0},B:1{1}):1{2},C:1{3}){4};',
      placements: [
        { p: [[1, 0.2, 0.1, 0.01], [0, 0.8, 0.5, 0.02]], n: ['q1'] },
        { p: [[3, 1, 0.3, 0.1]], nm: [['q2', 3], ['q3', 2]] },
      ],
    }
    const t = parseJplace(JSON.stringify(doc))
    const ix = indexTree(t)
    expect(ix.edgeByNum.size).toBe(4)
    expect(t.pqueries[0].places[0].edgeNum).toBe(0)
    expect(t.pqueries[0].mult).toBe(1)
    expect(t.pqueries[1].names).toEqual(['q2', 'q3'])
    expect(t.pqueries[1].mult).toBe(5)
  })
})

describe('split identity', () => {
  const src = parseNewick('((A:1,B:2)90:0.5,(C:1,(D:1,E:4)70:1)80:0.2,F:3);')
  const ix = indexTree(src)

  test('every branch has a distinct split id', () => {
    expect(ix.splitEdge.size).toBe(src.nodes.length - 1)
  })

  test('rerooting keeps the split set and total length', () => {
    const base = buildDisplay(ix, foldOps([]), ix.supportBySplit)
    const op = midpointOp(ix)!
    const re = buildDisplay(ix, foldOps([op]), ix.supportBySplit)
    const splits = (t: typeof base) => new Set(t.nodes.filter(d => d.split).map(d => d.split))
    expect(splits(re)).toEqual(splits(base))
    const total = (t: typeof base) => t.nodes.reduce((a, d) => a + d.len, 0)
    expect(total(re)).toBeCloseTo(total(base), 10)
  })

  test('Newick round-trip keeps splits and support', () => {
    const op = midpointOp(ix)!
    const d = buildDisplay(ix, foldOps([op]), ix.supportBySplit)
    const back = parseNewick(toNewick(d, 'support', true))
    expect(splitSet(back)).toEqual(new Set(ix.splitEdge.keys()))
    expect(new Map(indexTree(back).supportBySplit)).toEqual(new Map(ix.supportBySplit))
  })

  test('names and support together survive the NHX export', () => {
    const split = ix.split[src.nodes.findIndex(n => n.support === '70')]
    const st = foldOps([{ op: 'rename', key: `split:${split}`, name: 'Inner' }])
    const d = buildDisplay(ix, st, ix.supportBySplit)
    const back = parseNewick(toNewick(d, 'both', true))
    const inner = back.nodes.find(n => n.name === 'Inner')!
    expect(inner.support).toBe('70')
  })
})

describe('midpointOp', () => {
  test('roots at the middle of the longest tip-to-tip path', () => {
    // Longest path is A or B to D: 1 + 1 + 1 + 5 = 8, so the midpoint is 4 from D,
    // on D's own branch. D's side is the canonical side (it lacks the reference tip A).
    const src = parseNewick('((A:1,B:1):1,(C:1,D:5):1);')
    const ix = indexTree(src)
    const op = midpointOp(ix)!
    expect(op.op).toBe('reroot')
    if (op.op !== 'reroot') return
    expect(op.split).toBe(ix.split[byName(src, 'D')])
    expect(op.at).toBeCloseTo(0.8, 10)
  })
})

describe('op log', () => {
  const src = parseNewick('((A:1,B:1):1,(C:1,(D:1,E:1):1):1,F:1);')
  const ix = indexTree(src)
  const d = buildDisplay(ix, foldOps([]), new Map())
  const ab = d.nodes.find(n => n.tips === 2 && n.children.some(c => c.name === 'A'))!
  const de = d.nodes.find(n => n.tips === 2 && n.children.some(c => c.name === 'D'))!
  const ops: TreeOp[] = [
    { op: 'collapse', clade: ab.clade },
    { op: 'collapse', clade: de.clade },
    { op: 'expand', clade: de.clade },
    { op: 'rotate', clade: de.clade },
    { op: 'rename', key: 'tip:A', name: 'Alpha' },
    { op: 'rename', key: 'tip:A', name: 'Alpha2' },
    { op: 'rename', key: `split:${ab.split}`, name: 'AB', side: ab.clade },
    { op: 'style', clade: ab.clade, style: { bold: true, inherit: true } },
    { op: 'reroot', split: ab.split, at: 0.5 },
  ]

  test('stateOps reproduces the folded state', () => {
    const st = foldOps(ops)
    expect(stateJson(foldOps(stateOps(st)))).toBe(stateJson(st))
    expect(stateOps(st).length).toBeLessThan(ops.length)
  })

  test('a batch folds like its ops in sequence', () => {
    expect(stateJson(foldOps([{ op: 'batch', from: 'x', ops }]))).toBe(stateJson(foldOps(ops)))
  })

  test('matchOps keeps everything on a tree with the same tips', () => {
    const other = indexTree(parseNewick('((A:3,B:3):2,(C:1,(D:2,E:2):1):1,F:9);'))
    const { matched, unmatched } = matchOps(other, ops)
    expect(unmatched).toEqual([])
    expect(stateJson(foldOps(matched))).toBe(stateJson(foldOps(stateOps(foldOps(ops)))))
  })

  test('matchOps sets aside edits that refer to missing tips or branches', () => {
    const other = indexTree(parseNewick('((Z:1,B:1):1,(C:1,(D:1,E:1):1):1,F:1);'))
    const { unmatched } = matchOps(other, ops)
    expect(unmatched.some(o => o.op === 'rename' && o.key === 'tip:A')).toBe(true)
    expect(unmatched.some(o => o.op === 'collapse')).toBe(true)
  })
})

describe('single-child splicing', () => {
  // A chain of two unary nodes above AB: u2 (length 3) over u1 (2) over ab (1).
  // All three branches are drawn as one of length 6, and a point on each maps to
  // its own offset from the top of that merged branch.
  const src = parseNewick('((((A:1,B:1)ab:1)u1:2)u2:3,C:1,D:1);')
  const ix = indexTree(src)
  const d = buildDisplay(ix, foldOps([]), new Map())
  const on = (name: string, distal: number) => {
    const edge = byName(src, name)
    return pointOnBranch(d.byEdge.get(edge)![0], distal)
  }

  test('the merged branch keeps the full length', () => {
    const ab = d.nodes.find(n => n.tips === 2)!
    expect(ab.len).toBeCloseTo(6, 10)
    expect(ab.above).toBeCloseTo(5, 10)
  })

  test('points on each spliced segment land at their own offset', () => {
    expect(on('u2', 1).fromParent).toBeCloseTo(2, 10)
    expect(on('u1', 0.5).fromParent).toBeCloseTo(4.5, 10)
    expect(on('ab', 0.25).fromParent).toBeCloseTo(5.75, 10)
    expect(on('u2', 1).node).toBe(on('ab', 0.25).node)
  })
})
