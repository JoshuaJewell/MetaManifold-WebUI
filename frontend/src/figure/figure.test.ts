// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { describe, expect, test } from 'bun:test'
import { encodeTiff, lzwEncode, pngWithDpi } from './encode'
import { emptyHistory, record, redo, undo, undoKey } from './history'
import { legendEntries, mergeEntries, paneLegendEntries } from './legend'
import { PX_PER_PT, captionSplit, emptyPane, figureGeometry, textFrame, groupLetters, layoutText, newFigure, newGroup, parseRuns, resizeGroup, type FigureDoc } from './types'

const doc = (patch: Partial<FigureDoc> = {}): FigureDoc => ({ id: 'abcd1234', ...newFigure('Test'), ...patch })

describe('layout', () => {
  test('groups share the page height by weight, after margins and gaps', () => {
    const d = doc({ groups: [newGroup(1, 2), newGroup(2, 2), newGroup(1, 2)] })
    const g = figureGeometry(d)
    const heights = g.groups.map(x => x.box.h)
    expect(heights[1]).toBeCloseTo(2 * heights[0])
    expect(heights[0] + heights[1] + heights[2] + 2 * d.groupGapMm).toBeCloseTo(297 - 2 * d.page.marginMm)
    expect(g.groups[1].panes[3]).toEqual({
      x: 10 + (190 / 2), y: g.groups[1].box.y + heights[1] / 2, w: 190 / 2, h: heights[1] / 2,
    })
  })

  test('journal widths keep their width and take the chosen height', () => {
    const d = doc({ page: { size: 'single', widthMm: 999, heightMm: 70, marginMm: 0 } })
    expect(figureGeometry(d).page).toEqual({ w: 85, h: 70 })
  })

  test('letters skip groups with their own label', () => {
    const d = doc({ groups: [newGroup(), { ...newGroup(), label: 'S1' }, newGroup()] })
    expect(groupLetters(d)).toEqual(['a', 'S1', 'b'])
    expect(groupLetters({ ...d, letters: { ...d.letters, case: 'upper' } })).toEqual(['A', 'S1', 'B'])
  })

  test('letters can be turned off, keeping the ones a group sets', () => {
    const d = doc({ letters: { case: 'none', sizePt: 16 }, groups: [newGroup(), { ...newGroup(), label: 'S1' }] })
    expect(groupLetters(d)).toEqual(['', 'S1'])
  })

  test('column widths follow their weights', () => {
    const d = doc({ page: { size: 'custom', widthMm: 120, heightMm: 40, marginMm: 0 },
                    groups: [{ ...newGroup(1, 3), colWeights: [1, 2, 1] }] })
    const panes = figureGeometry(d).groups[0].panes
    expect(panes.map(p => p.w)).toEqual([30, 60, 30])
    expect(panes.map(p => p.x)).toEqual([0, 30, 90])
  })

  test('row heights follow their weights', () => {
    const d = doc({ page: { size: 'custom', widthMm: 100, heightMm: 60, marginMm: 0 },
                    groups: [{ ...newGroup(2, 1), rowWeights: [1, 2] }] })
    const panes = figureGeometry(d).groups[0].panes
    expect(panes.map(p => p.h)).toEqual([20, 40])
    expect(panes.map(p => p.y)).toEqual([0, 20])
  })

  test('resizing a grid keeps column weights and gives new columns 1', () => {
    const g = resizeGroup({ ...newGroup(1, 2), colWeights: [3, 2] }, 1, 3)
    expect(g.colWeights).toEqual([3, 2, 1])
    expect(resizeGroup(newGroup(1, 2), 1, 3).colWeights).toBeUndefined()
  })

  test('resizing a grid keeps panes at their row and column', () => {
    const g = newGroup(2, 2)
    g.panes[3] = { ...g.panes[3], title: 'bottom right' }
    const wider = resizeGroup(g, 2, 3)
    expect(wider.panes[4].title).toBe('bottom right')
    expect(wider.panes).toHaveLength(6)
    expect(resizeGroup(wider, 1, 1).panes).toHaveLength(1)
  })
})

describe('text panes', () => {
  test('markup splits into bold and italic runs', () => {
    expect(parseRuns('*Giardia* in **17**/35')).toEqual([
      { text: 'Giardia', bold: false, italic: true },
      { text: ' in ', bold: false, italic: false },
      { text: '17', bold: true, italic: false },
      { text: '/35', bold: false, italic: false },
    ])
    expect(parseRuns('***Hexamastix***')).toEqual([{ text: 'Hexamastix', bold: true, italic: true }])
    expect(parseRuns('2 * 3')).toEqual([{ text: '2 * 3', bold: false, italic: false }])
    expect(parseRuns('')).toEqual([{ text: '', bold: false, italic: false }])
  })

  test('a line too wide shrinks every line by the same factor', () => {
    const d = doc()
    // 10 px per character at any size: 'aaaaaaaaaa' is 100 px wide at its size, so it must shrink.
    const measure = (r: { text: string }, size: number) => r.text.length * size * 10 / (d.font.sizePt * PX_PER_PT)
    const lines = layoutText({ lines: [{ text: 'aaaaaaaaaa' }, { text: 'a', sizePt: d.font.sizePt * 2 }], align: 'left', valign: 'top' },
                             d, 50, 100, measure)
    const fit = 50 * 0.94 / 100
    expect(lines[0].sizePx).toBeCloseTo(d.font.sizePt * PX_PER_PT * fit)
    expect(lines[1].sizePx).toBeCloseTo(2 * d.font.sizePt * PX_PER_PT * fit)
  })

  test('a caption takes its lines from the bottom of a picture pane, at most half', () => {
    const d = doc()
    const pane = { ...emptyPane(), image: { file: 'a.png', fit: 'contain' as const },
                   text: { lines: [{ text: 'a' }, { text: 'b', sizePt: 20 }], align: 'center' as const, valign: 'top' as const } }
    const natural = (d.font.sizePt + 20) * PX_PER_PT * 1.2
    const split = captionSplit(pane, d, 200)
    expect(split.captionH).toBeCloseTo(natural)
    expect(split.imageH).toBeCloseTo(200 - natural)
    expect(captionSplit(pane, d, 40).captionH).toBe(20)
    expect(captionSplit({ ...pane, text: null }, d, 200)).toEqual({ imageH: 200, captionH: 0 })
  })

    test('turned text is laid out in the swapped box', () => {
    const t = { lines: [], align: 'center' as const, valign: 'middle' as const }
    expect(textFrame({ ...t, rotate: -90 }, 30, 100)).toEqual({ w: 100, h: 30, transform: 'translate(0,100) rotate(-90)' })
    expect(textFrame(t, 30, 100)).toEqual({ w: 30, h: 100, transform: undefined })
  })

  test('lines stack 1.2 em apart, centred in the pane', () => {
    const d = doc()
    const lines = layoutText({ lines: [{ text: 'a', sizePt: 20 }, { text: 'b' }], align: 'center', valign: 'middle' }, d, 200, 100)
    const big = 20 * PX_PER_PT, small = d.font.sizePt * PX_PER_PT
    const top = (100 - 1.2 * (big + small)) / 2
    expect(lines[0].y).toBeCloseTo(top + big * 0.95)
    expect(lines[1].y).toBeCloseTo(top + big * 1.2 + small * 0.95)
    expect(lines.map(l => [l.x, l.anchor])).toEqual([[100, 'middle'], [100, 'middle']])
  })
})

describe('legend', () => {
  const bar = (name: string, y: number[], extra = {}) => ({ type: 'bar', name, y, marker: { color: '#123456' }, ...extra })
  test('a bar series that never reaches 0.5% of the tallest stack is left out', () => {
    const fig = { data: [bar('Big', [0.9, 0.5]), bar('Small', [0.1, 0.004]), bar('Trace', [0.004, 0.001])], layout: {} }
    expect(mergeEntries([legendEntries(fig)]).map(e => e.name)).toEqual(['Big', 'Small'])
    // A name faint in the first chart keeps its place when a later chart shows it.
    const later = { data: [bar('Trace', [0.9]), bar('Big', [0.1])], layout: {} }
    expect(mergeEntries([legendEntries(fig), legendEntries(later)]).map(e => e.name)).toEqual(['Big', 'Small', 'Trace'])
    const turned = { data: [{ ...bar('Big', []), x: [1], orientation: 'h' }, { ...bar('Trace', []), x: [0.001], orientation: 'h' }], layout: {} }
    expect(mergeEntries([legendEntries(turned)]).map(e => e.name)).toEqual(['Big'])
  })

  test('genus names are italic, Other is not, and a pane can leave the legend', () => {
    const fig = { data: [bar('Hexamastix', [0.9]), bar('Other', [0.1])], layout: {} }
    const chart = { source: 'DADA2' as const, runs: [], table: 't', aggregate: false, kind: 'composition' as const, tag: 'rank' as const,
                    value: 'Genus', topN: 3, relative: true, mode: 'stacked' as const, keepEmpty: false, facet: null }
    const pane = { chart, title: null, xTitle: null }
    expect(paneLegendEntries(fig, pane).map(e => [e.name, e.italic ?? false])).toEqual([['Hexamastix', true], ['Other', false]])
    expect(paneLegendEntries(fig, { ...pane, style: { hideFromLegend: true } })).toEqual([])
  })
})

describe('encoders', () => {
  test('pngWithDpi inserts a valid pHYs chunk after IHDR', () => {
    const png = new Uint8Array(45)
    png.set([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a], 0)
    const out = pngWithDpi(png, 300)
    const v = new DataView(out.buffer)
    expect(String.fromCharCode(...out.subarray(37, 41))).toBe('pHYs')
    expect(v.getUint32(41)).toBe(Math.round(300 / 0.0254))
    expect(out[49]).toBe(1)
    expect(out.length).toBe(png.length + 21)
  })

  test('LZW output starts with a clear code and ends with EOI', () => {
    const bytes = lzwEncode(new Uint8Array([1, 1, 1, 1, 2, 2, 2, 2]))
    expect(bytes[0]).toBe(0x80)   // 256 in 9 bits, MSB first: 1000 0000 0...
    expect(bytes.length).toBeGreaterThan(2)
  })

  test('TIFF header and resolution', () => {
    const w = 3, h = 2
    const tif = encodeTiff(new Uint8ClampedArray(w * h * 4).fill(200), w, h, 600)
    const v = new DataView(tif.buffer)
    expect(String.fromCharCode(tif[0], tif[1])).toBe('II')
    expect(v.getUint16(2, true)).toBe(42)
    const ifd = v.getUint32(4, true)
    const tags = new Map<number, number>()
    for (let i = 0; i < v.getUint16(ifd, true); i++) {
      const p = ifd + 2 + i * 12
      tags.set(v.getUint16(p, true), v.getUint16(p + 2, true) === 3 ? v.getUint16(p + 8, true) : v.getUint32(p + 8, true))
    }
    expect(tags.get(256)).toBe(w)
    expect(tags.get(257)).toBe(h)
    expect(tags.get(259)).toBe(5)
    expect(v.getUint32(tags.get(282)!, true)).toBe(600)
  })
})

describe('undo', () => {
  test('edits undo and redo in order, and a new edit drops the redo steps', () => {
    let h = emptyHistory<string>()
    h = record(h, 'a', 1000)            // a -> b
    h = record(h, 'b', 2000)            // b -> c
    const u1 = undo(h, 'c')!
    expect(u1.doc).toBe('b')
    const u2 = undo(u1.history, 'b')!
    expect(u2.doc).toBe('a')
    expect(undo(u2.history, 'a')).toBeNull()
    const r1 = redo(u2.history, 'a')!
    expect(r1.doc).toBe('b')
    expect(record(r1.history, 'b', 5000).future).toEqual([])
  })

  test('edits close together undo as one step', () => {
    let h = emptyHistory<string>()
    h = record(h, 'a', 1000)            // a -> ab
    h = record(h, 'ab', 1200)           // ab -> abc, merged
    expect(h.past).toEqual(['a'])
    expect(undo(h, 'abc')!.doc).toBe('a')
  })

  test('keys', () => {
    const k = (key: string, mods: Partial<KeyboardEvent> = {}) =>
      undoKey({ key, ctrlKey: false, metaKey: false, shiftKey: false, altKey: false, ...mods })
    expect(k('z', { ctrlKey: true })).toBe('undo')
    expect(k('Z', { ctrlKey: true, shiftKey: true })).toBe('redo')
    expect(k('y', { ctrlKey: true })).toBe('redo')
    expect(k('z', { metaKey: true })).toBe('undo')
    expect(k('z')).toBeNull()
    expect(k('z', { ctrlKey: true, altKey: true })).toBeNull()
  })
})
