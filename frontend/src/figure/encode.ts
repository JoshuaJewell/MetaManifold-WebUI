// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).

const CRC_TABLE = (() => {
  const t = new Uint32Array(256)
  for (let n = 0; n < 256; n++) {
    let c = n
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1
    t[n] = c >>> 0
  }
  return t
})()

function crc32(bytes: Uint8Array): number {
  let c = 0xffffffff
  for (const b of bytes) c = CRC_TABLE[(c ^ b) & 0xff] ^ (c >>> 8)
  return (c ^ 0xffffffff) >>> 0
}

/** `png` with a pHYs chunk after IHDR, so the file carries its DPI. */
export function pngWithDpi(png: Uint8Array, dpi: number): Uint8Array<ArrayBuffer> {
  const ppm = Math.round(dpi / 0.0254)
  const chunk = new Uint8Array(21)
  const v = new DataView(chunk.buffer)
  v.setUint32(0, 9)
  chunk.set([0x70, 0x48, 0x59, 0x73], 4)
  v.setUint32(8, ppm)
  v.setUint32(12, ppm)
  chunk[16] = 1
  v.setUint32(17, crc32(chunk.subarray(4, 17)))
  // IHDR is the 8-byte signature plus a 25-byte chunk.
  const out = new Uint8Array(png.length + chunk.length)
  out.set(png.subarray(0, 33), 0)
  out.set(chunk, 33)
  out.set(png.subarray(33), 33 + chunk.length)
  return out
}

/** LZW, TIFF flavour: MSB-first codes, 9 to 12 bits, widening one code early. */
export function lzwEncode(input: Uint8Array): Uint8Array {
  const out: number[] = []
  let acc = 0, nbits = 0
  const put = (code: number, width: number) => {
    acc = (acc << width) | code
    nbits += width
    while (nbits >= 8) { out.push((acc >>> (nbits - 8)) & 0xff); nbits -= 8 }
    acc &= (1 << nbits) - 1
  }
  const CLEAR = 256, EOI = 257
  let dict = new Map<number, number>()
  let width = 9, next = 258
  put(CLEAR, width)
  if (input.length === 0) { put(EOI, width); if (nbits) out.push((acc << (8 - nbits)) & 0xff); return new Uint8Array(out) }
  let w = input[0]
  for (let i = 1; i < input.length; i++) {
    const k = input[i]
    const key = (w << 8) | k
    const hit = dict.get(key)
    if (hit !== undefined) { w = hit; continue }
    put(w, width)
    dict.set(key, next++)
    if (next === (1 << width) && width < 12) width++
    if (next === 4094) {
      put(CLEAR, width)
      dict = new Map()
      width = 9
      next = 258
    }
    w = k
  }
  put(w, width)
  next++
  if (next === (1 << width) && width < 12) width++
  put(EOI, width)
  if (nbits) out.push((acc << (8 - nbits)) & 0xff)
  return new Uint8Array(out)
}

/** Baseline RGB TIFF, LZW with horizontal differencing, resolution in DPI. */
export function encodeTiff(rgba: Uint8ClampedArray, width: number, height: number, dpi: number): Uint8Array<ArrayBuffer> {
  const rowBytes = width * 3
  const rowsPerStrip = Math.max(1, Math.floor(65536 / rowBytes))
  const strips: Uint8Array[] = []
  for (let y0 = 0; y0 < height; y0 += rowsPerStrip) {
    const rows = Math.min(rowsPerStrip, height - y0)
    const raw = new Uint8Array(rows * rowBytes)
    for (let r = 0; r < rows; r++) {
      const src = (y0 + r) * width * 4, dst = r * rowBytes
      for (let x = 0; x < width; x++) {
        raw[dst + x * 3]     = rgba[src + x * 4]
        raw[dst + x * 3 + 1] = rgba[src + x * 4 + 1]
        raw[dst + x * 3 + 2] = rgba[src + x * 4 + 2]
      }
      for (let i = rowBytes - 1; i >= 3; i--) raw[dst + i] = (raw[dst + i] - raw[dst + i - 3]) & 0xff
    }
    strips.push(lzwEncode(raw))
  }

  const entries: [number, number, number, number | number[]][] = []  // tag, type, count, value
  const SHORT = 3, LONG = 4, RATIONAL = 5
  const nEntries = 14
  const ifdOffset = 8
  const ifdSize = 2 + nEntries * 12 + 4
  let extra = ifdOffset + ifdSize
  const bpsOffset = extra; extra += 6
  const xresOffset = extra; extra += 8
  const yresOffset = extra; extra += 8
  const offsetsOffset = extra; extra += strips.length > 1 ? strips.length * 4 : 0
  const countsOffset = extra; extra += strips.length > 1 ? strips.length * 4 : 0
  let dataOffset = extra
  const stripOffsets = strips.map(s => { const o = dataOffset; dataOffset += s.length; return o })

  entries.push([256, LONG, 1, width], [257, LONG, 1, height], [258, SHORT, 3, bpsOffset],
               [259, SHORT, 1, 5], [262, SHORT, 1, 2],
               [273, LONG, strips.length, strips.length > 1 ? offsetsOffset : stripOffsets[0]],
               [277, SHORT, 1, 3], [278, LONG, 1, rowsPerStrip],
               [279, LONG, strips.length, strips.length > 1 ? countsOffset : strips[0].length],
               [282, RATIONAL, 1, xresOffset], [283, RATIONAL, 1, yresOffset],
               [284, SHORT, 1, 1], [296, SHORT, 1, 2], [317, SHORT, 1, 2])

  const buf = new Uint8Array(dataOffset)
  const v = new DataView(buf.buffer)
  buf.set([0x49, 0x49, 42, 0])
  v.setUint32(4, ifdOffset, true)
  v.setUint16(ifdOffset, nEntries, true)
  entries.forEach(([tag, type, count, value], i) => {
    const p = ifdOffset + 2 + i * 12
    v.setUint16(p, tag, true)
    v.setUint16(p + 2, type, true)
    v.setUint32(p + 4, count, true)
    if (type === SHORT && count === 1) v.setUint16(p + 8, value as number, true)
    else v.setUint32(p + 8, value as number, true)
  })
  v.setUint32(ifdOffset + 2 + nEntries * 12, 0, true)
  for (let i = 0; i < 3; i++) v.setUint16(bpsOffset + i * 2, 8, true)
  for (const o of [xresOffset, yresOffset]) { v.setUint32(o, Math.round(dpi), true); v.setUint32(o + 4, 1, true) }
  if (strips.length > 1) strips.forEach((s, i) => {
    v.setUint32(offsetsOffset + i * 4, stripOffsets[i], true)
    v.setUint32(countsOffset + i * 4, s.length, true)
  })
  strips.forEach((s, i) => buf.set(s, stripOffsets[i]))
  return buf
}
