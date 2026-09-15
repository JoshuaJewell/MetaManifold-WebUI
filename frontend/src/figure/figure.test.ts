// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { describe, expect, test } from 'bun:test'
import { encodeTiff, lzwEncode, pngWithDpi } from './encode'
import { figureGeometry, groupLetters, newFigure, newGroup, resizeGroup, type FigureDoc } from './types'

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

  test('resizing a grid keeps panes at their row and column', () => {
    const g = newGroup(2, 2)
    g.panes[3] = { ...g.panes[3], title: 'bottom right' }
    const wider = resizeGroup(g, 2, 3)
    expect(wider.panes[4].title).toBe('bottom right')
    expect(wider.panes).toHaveLength(6)
    expect(resizeGroup(wider, 1, 1).panes).toHaveLength(1)
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
