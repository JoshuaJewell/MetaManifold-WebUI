// SPDX-License-Identifier: MPL-2.0
import { describe, expect, test } from 'bun:test'
import { bestLwr, buildDisplay, foldOps, indexTree, parseJplace } from './model'
import { DEFAULT_SETTINGS, layoutTree } from './layout'
import { mergeSettings, readDoc } from './viewDoc'
import { fmtLwr } from './ui'

const tree = '((A:1{0},B:1{1}):1{2},C:1{3},D:1{4}){5};'
const source = (p: unknown[][], fields = ['edge_num', 'like_weight_ratio']) => parseJplace(JSON.stringify({
  tree, fields, placements: [{ n: ['q'], p }],
}))
const layout = (src: ReturnType<typeof source>, placements: 'best' | 'all') => {
  const ix = indexTree(src)
  return layoutTree(buildDisplay(ix, foldOps([]), new Map()), ix, src.pqueries, { ...DEFAULT_SETTINGS, placements })
}

describe('best placement LWR', () => {
  test('uses the maximum reported weight without rounding or normalizing', () => {
    const src = source([[0, 0.4], [1, 0.6]])
    expect(bestLwr(src.pqueries[0])).toBe(0.6)
    expect(src.pqueries[0].places[0].edgeNum).toBe(1)
    expect(bestLwr(source([[0, 0]]).pqueries[0])).toBe(0)
    expect(fmtLwr(0)).toBe('0')
  })

  test('absent field and absent candidate weights stay missing', () => {
    const noField = source([[0]], ['edge_num'])
    expect(noField.pqueries[0].places[0].lwr).toBeNull()
    expect(bestLwr(noField.pqueries[0])).toBeNull()
    expect(bestLwr(source([[0]]).pqueries[0])).toBeNull()
    expect(bestLwr(source([]).pqueries[0])).toBeNull()
    expect(fmtLwr(null)).toBe('Missing')
  })

  test('invalid values cannot become an apparent best weight', () => {
    for (const value of [null, '1', true, -0.1, 1.1]) {
      const q = source([[0, value], [1, 0.6]]).pqueries[0]
      expect(q.places.some(p => p.lwr === null)).toBe(true)
      expect(bestLwr(q)).toBeNull()
    }
  })

  test('weights follow field order and preserve tied candidates', () => {
    const src = source([[0.5, 0], [0.5, 1]], ['like_weight_ratio', 'edge_num'])
    expect(bestLwr(src.pqueries[0])).toBe(0.5)
    expect(src.pqueries[0].places).toHaveLength(2)
  })

  test('missing weights are not replaced by weighted placement counts', () => {
    const src = source([[0, 0.6], [1, null]])
    expect(layout(src, 'all').placements).toHaveLength(2)
    expect(layout(src, 'all').symbols).toHaveLength(0)
    expect(layout(src, 'best').placements).toHaveLength(0)
    const single = layout(source([[0]], ['edge_num']), 'best')
    expect(single.placements[0].lwr).toBeNull()
    expect(single.symbols[0].value).toBe(1) // one query, not an invented LWR
    expect(Number.isFinite(single.symbols[0].x)).toBe(true)
  })

  test('view changes do not change the raw best LWR', () => {
    const src = source([[0, 0.4], [1, 0.6]])
    const ix = indexTree(src)
    buildDisplay(ix, foldOps([{ op: 'reroot', split: ix.split[1], at: 0.5 }]), new Map())
    expect(bestLwr(src.pqueries[0])).toBe(0.6)
  })
})

describe('best LWR threshold setting', () => {
  test('older views get the default and saved thresholds survive loading', () => {
    expect(mergeSettings({}).bestLwrThreshold).toBe(0.95)
    expect(readDoc({ settings: { bestLwrThreshold: 0.8 } }).settings!.bestLwrThreshold).toBe(0.8)
  })

  test('invalid saved settings cannot disable comparison silently', () => {
    for (const bestLwrThreshold of [null, '0.5', NaN, Infinity]) {
      expect(mergeSettings({ bestLwrThreshold }).bestLwrThreshold).toBe(0.95)
    }
    expect(mergeSettings({ bestLwrThreshold: -1 }).bestLwrThreshold).toBe(0)
    expect(mergeSettings({ bestLwrThreshold: 2 }).bestLwrThreshold).toBe(1)
  })
})
