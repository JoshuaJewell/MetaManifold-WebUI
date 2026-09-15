// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import Plotly from 'plotly.js-dist-min'
import { PX_PER_MM, PX_PER_PT, MM_PER_IN, figureGeometry, titlePt, type FigureDoc } from './types'
import { layoutLegend, legendEntries, mergeEntries } from './legend'
import type { PlotSpec } from './charts'
import { encodeTiff, pngWithDpi } from './encode'

const SVG_NS = 'http://www.w3.org/2000/svg'

/** Plotly's SVG of one pane. */
async function paneSvg(fig: PlotSpec): Promise<SVGSVGElement> {
  const w = fig.layout.width as number, h = fig.layout.height as number
  const url = await Plotly.toImage({ data: fig.data as Plotly.Data[], layout: fig.layout as Partial<Plotly.Layout> },
                                   { format: 'svg', width: w, height: h })
  const text = decodeURIComponent(url.slice(url.indexOf(',') + 1))
  return new DOMParser().parseFromString(text, 'image/svg+xml').documentElement as unknown as SVGSVGElement
}

// Pane SVGs reuse ids such as clip paths, so each pane's ids get their own prefix.
function prefixIds(root: Element, prefix: string) {
  const all = [root, ...Array.from(root.querySelectorAll('*'))]
  for (const el of all) if (el.id) el.id = prefix + el.id
  const ref = (v: string) => v
    .replace(/url\((['"]?)#([^'")]+)\1\)/g, (_m, q, id) => `url(${q}#${prefix}${id}${q})`)
  for (const el of all) {
    for (const attr of Array.from(el.attributes)) {
      if ((attr.name === 'href' || attr.name === 'xlink:href') && attr.value.startsWith('#')) {
        el.setAttribute(attr.name, '#' + prefix + attr.value.slice(1))
      } else if (attr.value.includes('url(')) {
        el.setAttribute(attr.name, ref(attr.value))
      }
    }
  }
}

/** The whole page as one SVG, in CSS pixels at 96 per inch. */
export async function figureSvg(doc: FigureDoc, panes: (PlotSpec | null)[][]): Promise<SVGSVGElement> {
  const geo = figureGeometry(doc)
  const W = geo.page.w * PX_PER_MM, H = geo.page.h * PX_PER_MM
  const svg = document.createElementNS(SVG_NS, 'svg')
  svg.setAttribute('xmlns', SVG_NS)
  svg.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink')
  svg.setAttribute('width', `${geo.page.w}mm`)
  svg.setAttribute('height', `${geo.page.h}mm`)
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`)
  const rect = (x: number, y: number, w: number, h: number, fill: string) => {
    const r = document.createElementNS(SVG_NS, 'rect')
    for (const [k, v] of Object.entries({ x, y, width: w, height: h, fill })) r.setAttribute(k, String(v))
    return r
  }
  svg.appendChild(rect(0, 0, W, H, '#ffffff'))

  const family = doc.font.family === 'Arial' ? 'Arial, Helvetica, sans-serif' : '"Times New Roman", Times, serif'
  let n = 0
  for (const [gi, g] of geo.groups.entries()) {
    const group = doc.groups[gi]
    const b = g.box
    if (group.background) svg.appendChild(rect(b.x * PX_PER_MM, b.y * PX_PER_MM, b.w * PX_PER_MM, b.h * PX_PER_MM, group.background))
    for (const [pi, box] of g.panes.entries()) {
      const fig = panes[gi]?.[pi]
      if (!fig) continue
      const pane = await paneSvg(fig)
      prefixIds(pane, `f${n++}-`)
      const wrap = document.createElementNS(SVG_NS, 'g')
      wrap.setAttribute('transform', `translate(${box.x * PX_PER_MM},${box.y * PX_PER_MM})`)
      for (const child of Array.from(pane.childNodes)) wrap.appendChild(document.importNode(child, true))
      svg.appendChild(wrap)
    }
    const text = (x: number, y: number, size: number, content: string, anchor = 'start') => {
      const el = document.createElementNS(SVG_NS, 'text')
      for (const [k, v] of Object.entries({ x, y, 'font-family': family, 'font-size': size, fill: '#444444', 'text-anchor': anchor }))
        el.setAttribute(k, String(v))
      el.textContent = content
      svg.appendChild(el)
    }
    if (g.title && group.title) {
      const size = titlePt(doc) * PX_PER_PT
      text((g.title.x + g.title.w / 2) * PX_PER_MM, (g.title.y + g.title.h / 2) * PX_PER_MM + size * 0.35, size, group.title, 'middle')
    }
    if (g.legend) {
      const entries = mergeEntries((panes[gi] ?? []).filter((p): p is PlotSpec => p != null).map(legendEntries))
      const lx = g.legend.x * PX_PER_MM, ly = g.legend.y * PX_PER_MM, lh = g.legend.h * PX_PER_MM
      const lay = layoutLegend(entries, g.legend.w * PX_PER_MM, doc.font.sizePt * PX_PER_PT, family)
      for (const it of lay.items) {
        const cy = ly + lh / 2
        const sw = document.createElementNS(SVG_NS, it.entry.shape === 'circle' ? 'circle' : 'rect')
        const attrs: Record<string, string | number> = it.entry.shape === 'circle'
          ? { cx: lx + it.x + lay.swatch / 2, cy, r: lay.swatch / 2 }
          : { x: lx + it.x, y: cy - lay.swatch / 2, width: lay.swatch, height: lay.swatch }
        attrs.fill = it.entry.fill
        if (it.entry.stroke) { attrs.stroke = it.entry.stroke; attrs['stroke-width'] = 1 }
        for (const [k, v] of Object.entries(attrs)) sw.setAttribute(k, String(v))
        svg.appendChild(sw)
        text(lx + it.textX, cy + lay.fontPx * 0.35, lay.fontPx, it.entry.name)
      }
    }
    const size = doc.letters.sizePt * PX_PER_PT
    const t = document.createElementNS(SVG_NS, 'text')
    t.setAttribute('x', String(b.x * PX_PER_MM + size * 0.15))
    t.setAttribute('y', String(b.y * PX_PER_MM + size * 0.9))
    t.setAttribute('font-family', family)
    t.setAttribute('font-size', String(size))
    t.setAttribute('font-weight', 'bold')
    t.setAttribute('fill', '#000000')
    t.textContent = g.letter
    svg.appendChild(t)
  }
  return svg
}

// The PDF's built-in fonts cover Windows-1252 only.
const PDF_TEXT: [RegExp, string][] = [[/\u2212/g, '-'], [/\u2264/g, '<='], [/\u2265/g, '>='], [/\u03bc/g, '\u00b5'], [/\u2009|\u202f/g, ' ']]

export async function figurePdf(doc: FigureDoc, svg: SVGSVGElement): Promise<Blob> {
  const [{ jsPDF }, { svg2pdf }] = await Promise.all([import('jspdf'), import('svg2pdf.js')])
  const geo = figureGeometry(doc)
  const clone = svg.cloneNode(true) as SVGSVGElement
  const walker = document.createTreeWalker(clone, NodeFilter.SHOW_TEXT)
  for (let node = walker.nextNode(); node; node = walker.nextNode())
    for (const [re, to] of PDF_TEXT) node.nodeValue = node.nodeValue!.replace(re, to)
  const pdf = new jsPDF({ unit: 'mm', format: [geo.page.w, geo.page.h],
                          orientation: geo.page.w > geo.page.h ? 'landscape' : 'portrait' })
  // svg2pdf reads computed styles, so the SVG has to be in the document while it runs.
  const host = document.createElement('div')
  host.style.cssText = 'position:fixed;left:-100000px;top:0;visibility:hidden'
  host.appendChild(clone)
  document.body.appendChild(host)
  // svg2pdf sizes em offsets from the wrong font, which shifts Plotly's titles.
  for (const el of Array.from(clone.querySelectorAll('text, tspan'))) {
    const size = parseFloat(getComputedStyle(el).fontSize)
    for (const attr of ['dx', 'dy']) {
      const v = el.getAttribute(attr)
      if (v?.endsWith('em')) el.setAttribute(attr, String(parseFloat(v) * size))
    }
  }
  try {
    await svg2pdf(clone, pdf, { x: 0, y: 0, width: geo.page.w, height: geo.page.h })
  } finally {
    host.remove()
  }
  pdf.setProperties({ title: doc.title })
  return pdf.output('blob')
}

/** The page drawn onto a white canvas at `dpi`. */
async function rasterise(doc: FigureDoc, svg: SVGSVGElement, dpi: number): Promise<HTMLCanvasElement> {
  const geo = figureGeometry(doc)
  const w = Math.round(geo.page.w / MM_PER_IN * dpi), h = Math.round(geo.page.h / MM_PER_IN * dpi)
  const clone = svg.cloneNode(true) as SVGSVGElement
  clone.setAttribute('width', String(w))
  clone.setAttribute('height', String(h))
  const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml' }))
  try {
    const img = new Image()
    img.decoding = 'sync'
    await new Promise<void>((ok, fail) => { img.onload = () => ok(); img.onerror = () => fail(new Error('Could not draw the figure')); img.src = url })
    const canvas = document.createElement('canvas')
    canvas.width = w
    canvas.height = h
    const ctx = canvas.getContext('2d')!
    ctx.fillStyle = '#ffffff'
    ctx.fillRect(0, 0, w, h)
    ctx.drawImage(img, 0, 0, w, h)
    return canvas
  } finally {
    URL.revokeObjectURL(url)
  }
}

export async function figurePng(doc: FigureDoc, svg: SVGSVGElement, dpi: number): Promise<Blob> {
  const canvas = await rasterise(doc, svg, dpi)
  const blob = await new Promise<Blob>((ok, fail) => canvas.toBlob(b => b ? ok(b) : fail(new Error('PNG encoding failed')), 'image/png'))
  return new Blob([pngWithDpi(new Uint8Array(await blob.arrayBuffer()), dpi)], { type: 'image/png' })
}

export async function figureTiff(doc: FigureDoc, svg: SVGSVGElement, dpi: number): Promise<Blob> {
  const canvas = await rasterise(doc, svg, dpi)
  const { width, height } = canvas
  const rgba = canvas.getContext('2d')!.getImageData(0, 0, width, height).data
  return new Blob([encodeTiff(rgba, width, height, dpi)], { type: 'image/tiff' })
}
