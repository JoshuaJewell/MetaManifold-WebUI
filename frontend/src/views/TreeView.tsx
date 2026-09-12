// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import { AddToReport } from '../components/AddToReport'
import { Link, useParams } from 'react-router-dom'
import { api } from '../api/client'
import type { TreeFile, TreeFileInfo } from '../api/types'
import { useApi } from '../hooks/useApi'
import { Skeleton } from '../components/Skeleton'
import { useToast } from '../components/Toast'
import {
  buildDisplay, foldOps, indexTree, mapSupport, midpointOp, nodeKey, parseTreeFile, pointOnBranch,
  rerootAbove, tipsBelow, toNewick, unmatchedOps,
  type DNode, type SourceTree, type TreeOp,
} from '../tree/model'
import { DEFAULT_SETTINGS, layoutTree, symbolPath, type Dataset, type SymbolStyle, type TreeSettings } from '../tree/layout'
import styles from '../tree/Tree.module.css'
import { mergeSettings, readDoc, type TreeViewDoc } from '../tree/viewDoc'
import { fmt, PLACE_SELECTED } from '../tree/ui'
import { TreeShapes } from '../tree/TreeShapes'
import { Legend, legendItems } from '../tree/Legend'
import { ImportPanel } from '../tree/ImportPanel'
import { DatasetControls } from '../tree/controls'
import { Inspector } from '../tree/Inspector'
import { saveBlob } from '../utils/download'

function renameOp(key: string, name: string | null, side: string): TreeOp {
  return key.startsWith('split:') && name !== null ? { op: 'rename', key, name, side } : { op: 'rename', key, name }
}

/** Where a tree and its saved view come from: a study's tree list or the reference tree library. */
export interface TreeSource {
  key: string
  back: { to: string; label: string }
  get: (file: string) => Promise<TreeFile>
  saveView: (file: string, doc: unknown) => Promise<unknown>
  /** Other trees support and edits can be imported from. */
  list: () => Promise<TreeFileInfo[]>
  /** The study whose report a figure goes to, if any. */
  study: string | null
}

export function TreeView() {
  const { study, file } = useParams<{ study: string; file: string }>()
  const source = useMemo<TreeSource>(() => ({
    key: study!,
    back: { to: `/${study}?view=trees`, label: 'All trees' },
    get: f => api.trees.get(study!, f),
    saveView: (f, doc) => api.trees.saveView(study!, f, doc),
    list: () => api.trees.list(study!),
    study: study!,
  }), [study])
  return <TreeLoader source={source} file={file!} />
}

/** A reference tree from the library, with its view saved beside it. */
export function ReferenceTreeView() {
  const { id } = useParams<{ id: string }>()
  const source = useMemo<TreeSource>(() => ({
    key: `ref:${id}`,
    back: { to: `/reference-trees?tree=${id}`, label: 'Reference tree' },
    get: () => api.referenceTrees.tree(id!),
    saveView: (_, doc) => api.referenceTrees.saveView(id!, doc),
    list: () => Promise.resolve([]),
    study: null,
  }), [id])
  return <TreeLoader source={source} file="" />
}

function TreeLoader({ source, file }: { source: TreeSource; file: string }) {
  const fetcher = useCallback(() => source.get(file), [source, file])
  const { data, loading, error } = useApi(fetcher)

  if (loading) return <Skeleton lines={6} />
  if (error || !data) return <p className="error-msg">{error ?? 'Tree not found'}</p>
  return <TreeEditor key={`${source.key}/${file}`} source={source} tree={data} />
}

function TreeEditor({ source, tree }: { source: TreeSource; tree: TreeFile }) {
  const toast = useToast()
  const parsed = useMemo(() => {
    try { return { src: parseTreeFile(tree.format, tree.content), error: null } }
    catch (e) { return { src: null, error: String((e as Error).message ?? e) } }
  }, [tree])

  if (!parsed.src) return (
    <>
      <Header source={source} file={tree.file} />
      <p className="error-msg">Could not read {tree.file}: {parsed.error}</p>
    </>
  )
  return <TreeCanvas source={source} tree={tree} src={parsed.src} toast={toast} />
}

function Header({ source, file, children }: { source: TreeSource; file: string; children?: React.ReactNode }) {
  return (
    <div className="page-header" style={{ display: 'flex', alignItems: 'flex-start', justifyContent: 'space-between', gap: 16, marginBottom: 12 }}>
      <div>
        <h1 style={{ wordBreak: 'break-all' }}>{file}</h1>
        {children}
      </div>
      <Link className="btn" to={source.back.to}>{source.back.label}</Link>
    </div>
  )
}

interface Editing { key: string; side: string; value: string; left: number; top: number; width: number }

function TreeCanvas({ source, tree, src, toast }: {
  source: TreeSource; tree: TreeFile; src: SourceTree; toast: ReturnType<typeof useToast>
}) {
  const initial = useMemo(() => readDoc(tree.view), [tree])
  const ix = useMemo(() => indexTree(src), [src])
  const [hist, setHist] = useState<{ ops: TreeOp[]; redo: TreeOp[] }>({ ops: initial.ops ?? [], redo: [] })
  const { ops, redo } = hist
  const [settings, setSettings] = useState<TreeSettings>(initial.settings ?? DEFAULT_SETTINGS)
  const [selClade, setSelClade] = useState<string | null>(initial.selected ?? null)
  const [selQuery, setSelQuery] = useState<number | null>(null)
  const [editing, setEditing] = useState<Editing | null>(null)
  const [search, setSearch] = useState('')
  const [matchIdx, setMatchIdx] = useState(0)
  const [saveState, setSaveState] = useState<'saved' | 'pending' | 'error'>('saved')
  const scrollRef = useRef<HTMLDivElement>(null)
  const svgRef = useRef<SVGSVGElement>(null)
  const scrollPos = useRef(initial.scroll ?? { left: 0, top: 0 })

  const set = <K extends keyof TreeSettings>(k: K, v: TreeSettings[K]) => setSettings(s => ({ ...s, [k]: v }))
  const setSym = (ds: Dataset, patch: Partial<SymbolStyle>) =>
    setSettings(s => ({ ...s, symbols: { ...s.symbols, [ds]: { ...s.symbols[ds], ...patch } } }))

  const sourceChanged = !!initial.source && initial.source.sha256 !== tree.sha256
  const unmatched = useMemo(() => sourceChanged ? unmatchedOps(ix, initial.ops ?? []) : 0, [sourceChanged, ix, initial])

  // Consecutive style edits to one clade (a colour picker drag) count as one edit.
  const push = useCallback((op: TreeOp) => setHist(h => {
    const last = h.ops[h.ops.length - 1]
    const same = op.op === 'style' && last?.op === 'style' && last.clade === op.clade
    return { ops: same ? [...h.ops.slice(0, -1), op] : [...h.ops, op], redo: [] }
  }), [])
  const [importOpen, setImportOpen] = useState(false)
  const listFetcher = useCallback(() => source.list(), [source])
  const { data: studyTrees } = useApi(listFetcher)

  const state = useMemo(() => foldOps(ops), [ops])
  // Imported support overrides values written in this file.
  const support = useMemo(() => {
    if (!state.support.size) return ix.supportBySplit
    return new Map([...ix.supportBySplit, ...state.support])
  }, [ix, state.support])

  // Views saved while support was a live link to another file become an import.
  useEffect(() => {
    const f = settings.supportTree
    if (!f) return
    let cancelled = false
    source.get(f)
      .then(t => {
        if (cancelled) return
        const { support: values } = mapSupport(ix, parseTreeFile(t.format, t.content))
        push({ op: 'support', from: f, values: [...values] })
      })
      .catch(e => { if (!cancelled) toast.error(`Support from ${f}: ${e.message ?? e}`) })
      .finally(() => { if (!cancelled) setSettings(s => ({ ...s, supportTree: null })) })
    return () => { cancelled = true }
  }, [settings.supportTree, source, ix, push, toast])

  const display = useMemo(() => buildDisplay(ix, state, support), [ix, state, support])
  const layout = useMemo(() => layoutTree(display, ix, src.pqueries, settings), [display, ix, src, settings])

  const legendRows = legendItems(layout, settings)
  const legendH = legendRows.length ? legendRows.reduce((a, r) => a + r.height, 0) + 12 : 0

  const selected: DNode | null = useMemo(() => {
    if (!selClade) return null
    return display.nodes.find(d => d.clade === selClade) ?? null
  }, [display, selClade])

  // Undo moves the last op onto the redo stack.
  const undo = useCallback(() => setHist(h => h.ops.length
    ? { ops: h.ops.slice(0, -1), redo: [...h.redo, h.ops[h.ops.length - 1]] } : h), [])
  const redoOne = useCallback(() => setHist(h => h.redo.length
    ? { ops: [...h.ops, h.redo[h.redo.length - 1]], redo: h.redo.slice(0, -1) } : h), [])

  // Autosave the view document.
  const firstSave = useRef(true)
  const saveTimer = useRef<number>()
  const saveNow = useCallback(() => {
    const doc: TreeViewDoc = {
      version: 1,
      source: { file: tree.file, sha256: tree.sha256 },
      ops, settings, scroll: scrollPos.current, selected: selClade,
    }
    setSaveState('pending')
    source.saveView(tree.file, doc)
      .then(() => setSaveState('saved'))
      .catch(() => setSaveState('error'))
  }, [ops, settings, selClade, source, tree])
  // A save still waiting on its timer is sent when the page unmounts.
  const savePending = useRef(false)
  const latestSave = useRef(saveNow)
  latestSave.current = saveNow
  const flush = useCallback(() => { savePending.current = false; latestSave.current() }, [])
  useEffect(() => {
    if (firstSave.current) { firstSave.current = false; return }
    setSaveState('pending')
    savePending.current = true
    window.clearTimeout(saveTimer.current)
    saveTimer.current = window.setTimeout(flush, 600)
    return () => window.clearTimeout(saveTimer.current)
  }, [saveNow, flush])

  // Restore the scroll position once the first layout is in the DOM.
  useEffect(() => {
    const el = scrollRef.current
    if (!el) return
    el.scrollLeft = scrollPos.current.left
    el.scrollTop = scrollPos.current.top
  }, [])
  // Ctrl+wheel zooms; React's wheel listener is passive and cannot cancel the page zoom.
  useEffect(() => {
    const el = scrollRef.current
    if (!el) return
    const onWheel = (e: WheelEvent) => {
      if (!e.ctrlKey) return
      e.preventDefault()
      setSettings(s => ({ ...s, zoom: Math.min(10, Math.max(0.1, +(s.zoom * (e.deltaY < 0 ? 1.1 : 1 / 1.1)).toFixed(3))) }))
    }
    el.addEventListener('wheel', onWheel, { passive: false })
    return () => el.removeEventListener('wheel', onWheel)
  }, [])
  const scrollTimer = useRef<number>()
  const onScroll = () => {
    const el = scrollRef.current
    if (!el) return
    scrollPos.current = { left: el.scrollLeft, top: el.scrollTop }
    savePending.current = true
    window.clearTimeout(scrollTimer.current)
    scrollTimer.current = window.setTimeout(flush, 1200)
  }
  useEffect(() => () => {
    window.clearTimeout(saveTimer.current)
    window.clearTimeout(scrollTimer.current)
    if (savePending.current) flush()
  }, [flush])

  const scrollToNode = useCallback((id: number) => {
    const n = layout.byId.get(id)
    const el = scrollRef.current
    if (!n || !el) return
    el.scrollTo({ left: Math.max(0, n.x - el.clientWidth / 2), top: Math.max(0, n.y + legendH - el.clientHeight / 2), behavior: 'smooth' })
  }, [layout, legendH])

  const toggleCollapse = useCallback((d: DNode) => {
    if (d.tipKey !== null || !d.parent) return
    push(d.collapsed ? { op: 'expand', clade: d.clade } : { op: 'collapse', clade: d.clade })
  }, [push])

  const reroot = useCallback((d: DNode) => {
    const op = rerootAbove(ix, d)
    if (op) push(op)
  }, [ix, push])

  const startRename = useCallback((d: DNode, target?: Element) => {
    const key = nodeKey(d)
    const el = scrollRef.current
    const n = layout.byId.get(d.id)
    if (!key || !el || !n) return
    const box = el.getBoundingClientRect()
    let left = n.x + 4, top = n.y + legendH - 10
    if (target) {
      const r = target.getBoundingClientRect()
      left = r.left - box.left + el.scrollLeft
      top = r.top - box.top + el.scrollTop - 2
    }
    setEditing({ key, side: d.clade, value: d.name, left, top, width: Math.max(160, d.name.length * settings.fontSize * 0.62 + 24) })
  }, [layout, legendH, settings.fontSize])

  const commitRename = () => {
    if (!editing) return
    const current = state.names.get(editing.key) ?? null
    const value = editing.value.trim()
    const original = editing.key.startsWith('tip:') ? editing.key.slice(4) : ''
    const next = value === '' || value === original ? null : value
    if (next !== current) push(renameOp(editing.key, next, editing.side))
    setEditing(null)
  }

  // Search over tip names, including tips inside collapsed clades.
  const matches = useMemo(() => {
    const q = search.trim().toLowerCase()
    if (!q) return [] as DNode[]
    return display.nodes.filter(d => d.tipKey !== null && d.name.toLowerCase().includes(q))
  }, [search, display])
  const matchSet = useMemo(() => new Set(matches.map(m => m.id)), [matches])
  const visibleOwner = useCallback((d: DNode): DNode => {
    let top: DNode = d
    for (let p = d.parent; p; p = p.parent) if (p.collapsed) top = p
    return top
  }, [])
  const jumpMatch = (step: number) => {
    if (!matches.length) return
    const k = ((matchIdx + step) % matches.length + matches.length) % matches.length
    setMatchIdx(k)
    const owner = visibleOwner(matches[k])
    setSelClade(owner.clade)
    scrollToNode(owner.id)
  }

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const t = e.target as HTMLElement
      if (t && (t.tagName === 'INPUT' || t.tagName === 'SELECT' || t.tagName === 'TEXTAREA')) return
      const mod = e.ctrlKey || e.metaKey
      if (mod && e.key.toLowerCase() === 'z' && !e.shiftKey) { e.preventDefault(); undo() }
      else if (mod && (e.key.toLowerCase() === 'y' || (e.key.toLowerCase() === 'z' && e.shiftKey))) { e.preventDefault(); redoOne() }
      else if (!mod && selected) {
        if (e.key === 'c') toggleCollapse(selected)
        else if (e.key === 'r') reroot(selected)
        else if (e.key === 'o' && selected.tipKey === null) push({ op: 'rotate', clade: selected.clade })
        else if (e.key === 'F2') { e.preventDefault(); startRename(selected) }
        else if (e.key === 'Escape') setSelClade(null)
        else if (e.key === 'ArrowUp' && selected.parent) { e.preventDefault(); setSelClade(selected.parent.clade) }
      }
    }
    window.addEventListener('keydown', onKey)
    return () => window.removeEventListener('keydown', onKey)
  }, [undo, redoOne, selected, toggleCollapse, reroot, startRename, push])

  // Export
  const baseName = tree.file.replace(/\.[^.]+$/, '')
  const svgMarkup = (): string => {
    const svg = svgRef.current!.cloneNode(true) as SVGSVGElement
    svg.querySelectorAll('[data-ui]').forEach(el => el.remove())
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
    svg.setAttribute('color', '#000000')
    svg.removeAttribute('class')
    svg.style.cssText = ''
    const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect')
    bg.setAttribute('width', '100%'); bg.setAttribute('height', '100%'); bg.setAttribute('fill', '#ffffff')
    svg.insertBefore(bg, svg.firstChild)
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg)
  }
  const exportSvg = () => saveBlob(new Blob([svgMarkup()], { type: 'image/svg+xml' }), `${baseName}.svg`)
  const exportPng = () => {
    const img = new Image()
    const url = URL.createObjectURL(new Blob([svgMarkup()], { type: 'image/svg+xml' }))
    img.onload = () => {
      const k = 3
      const canvas = document.createElement('canvas')
      canvas.width = layout.width * k
      canvas.height = (layout.height + legendH) * k
      const ctx = canvas.getContext('2d')!
      ctx.scale(k, k)
      ctx.drawImage(img, 0, 0)
      URL.revokeObjectURL(url)
      canvas.toBlob(b => b ? saveBlob(b, `${baseName}.png`) : toast.error('PNG export failed: the image is too large for the browser canvas'), 'image/png')
    }
    img.onerror = () => { URL.revokeObjectURL(url); toast.error('PNG export failed') }
    img.src = url
  }
  const exportNewick = () => saveBlob(
    new Blob([toNewick(display, settings.newickInternal, src.hasLengths) + '\n'], { type: 'text/plain' }),
    `${baseName}.nwk`,
  )

  const selectedPlacements = useMemo(() => {
    if (!selected) return []
    return layout.placements.filter(p => p.node.id === selected.id)
  }, [layout, selected])

  const describe = useCallback((d: DNode): string => {
    if (d.tipKey !== null) return d.name
    if (d.name) return d.name
    const tips = tipsBelow(d)
    if (!tips.length) return 'root'
    return tips.length === 1 ? tips[0].name : `${tips.length} tips: ${tips[0].name} … ${tips[tips.length - 1].name}`
  }, [])

  const queryRows = useMemo(() => src.pqueries.map((q, qi) => {
    const best = q.places[0]
    const edge = best ? ix.edgeByNum.get(best.edgeNum) : undefined
    const onEdge = edge !== undefined ? display.byEdge.get(edge)?.[0] : undefined
    const dn = onEdge && best ? pointOnBranch(onEdge, best.distal).node : undefined
    return { qi, name: q.names.join(', '), lwr: best?.lwr ?? 0, n: q.places.length, where: dn ? describe(dn) : '', dn }
  }), [src, ix, display, describe])

  const selectQuery = (qi: number) => {
    setSelQuery(qi === selQuery ? null : qi)
    const row = queryRows[qi]
    if (row?.dn && qi !== selQuery) {
      const owner = visibleOwner(row.dn)
      setSelClade(owner.clade)
      scrollToNode(owner.id)
    }
  }

  const sel = selected ? layout.byId.get(selected.id) : undefined
  const selSubtree = useMemo(() => {
    const out = new Set<number>()
    if (!sel) return out
    const stack = [sel]
    while (stack.length) { const n = stack.pop()!; out.add(n.d.id); stack.push(...n.kids) }
    return out
  }, [sel])

  return (
    <>
      <Header source={source} file={tree.file}>
        <p style={{ fontSize: '.85rem', color: 'var(--color-muted-fg)' }}>
          {ix.tipCount} tips
          {src.pqueries.length > 0 && <> · {src.pqueries.length} placed queries</>}
          {' · '}{ops.length} edit{ops.length === 1 ? '' : 's'}
          {' · '}
          <span className={saveState === 'error' ? styles.saveError : undefined}>
            {saveState === 'saved' ? 'View saved' : saveState === 'pending' ? 'Saving…' : 'View not saved'}
          </span>
        </p>
      </Header>

      {importOpen && (
        <ImportPanel load={source.get} file={tree.file} ix={ix}
          candidates={(studyTrees ?? []).filter(t => t.file !== tree.file).map(t => t.file)}
          onClose={() => setImportOpen(false)}
          onApply={(from, copied, fromSettings) => {
            if (copied.length) push({ op: 'batch', from, ops: copied })
            if (fromSettings) setSettings(mergeSettings(fromSettings))
            setImportOpen(false)
            toast.success(`Imported from ${from}`)
          }} />
      )}

      {sourceChanged && (
        <p className={styles.notice}>
          {tree.file} has changed since this view was saved. Edits were replayed by branch identity
          {unmatched > 0 ? `; ${unmatched} refer to branches or tips no longer in the tree.` : ' and all of them still match.'}
        </p>
      )}

      <div className={`card ${styles.toolbar}`}>
        <div className={styles.row}>
          <span className={styles.label}>Layout</span>
          <select value={settings.layout} onChange={e => set('layout', e.target.value as TreeSettings['layout'])}>
            <option value="rectangular">Rectangular</option>
            <option value="slanted">Slanted</option>
            <option value="circular">Circular</option>
            <option value="unrooted">Unrooted</option>
          </select>
          {settings.layout === 'circular' && <>
            <label>Arc <input type="number" min={10} max={360} step={10} value={settings.arc} className={styles.num}
              onChange={e => set('arc', Math.min(360, Math.max(10, Number(e.target.value) || 360)))} />°</label>
            <label>Start <input type="number" min={-360} max={360} step={15} value={settings.arcStart} className={styles.num}
              onChange={e => set('arcStart', Number(e.target.value) || 0)} />°</label>
          </>}
          <label><input type="checkbox" checked={settings.branchLengths} disabled={!src.hasLengths}
            onChange={e => set('branchLengths', e.target.checked)} /> Branch lengths</label>
          {settings.branchLengths && src.hasLengths && <>
            <label>drawn <select value={settings.branchTransform}
              onChange={e => set('branchTransform', e.target.value as TreeSettings['branchTransform'])}>
              <option value="none">as in file</option>
              <option value="cap">capped</option>
              <option value="sqrt">square root</option>
              <option value="equal">equal</option>
            </select></label>
            {settings.branchTransform === 'cap' && (
              <label>at <input type="number" min={0} step="any" className={styles.num}
                placeholder={layout.capValue !== null ? fmt(layout.capValue, 3) : ''}
                value={settings.capAt ?? ''}
                onChange={e => set('capAt', e.target.value === '' ? null : Math.max(0, Number(e.target.value)))} /></label>
            )}
            {(settings.branchTransform === 'none' || settings.branchTransform === 'cap') && <>
              <label><input type="checkbox" checked={settings.scaleBar} onChange={e => set('scaleBar', e.target.checked)} /> Scale bar</label>
              {settings.scaleBar && (
                <input type="number" min={0} step="any" className={styles.num} title="Scale bar length" aria-label="Scale bar length"
                  placeholder={layout.scaleBar ? fmt(layout.scaleBar.value, 3) : ''}
                  value={settings.scaleBarLength ?? ''}
                  onChange={e => set('scaleBarLength', e.target.value === '' ? null : Math.max(0, Number(e.target.value)))} />
              )}
            </>}
          </>}
          <label>Ladderize <select value={settings.ladderize} onChange={e => set('ladderize', e.target.value as TreeSettings['ladderize'])}>
            <option value="none">Off</option>
            <option value="up">Small clades first</option>
            <option value="down">Large clades first</option>
          </select></label>
          {settings.ladderize !== 'none' && (
            <label>by <select value={settings.ladderizeBy} onChange={e => set('ladderizeBy', e.target.value as TreeSettings['ladderizeBy'])}>
              <option value="tips">tips</option>
              <option value="rows">visible rows</option>
            </select></label>
          )}
          <label>Font <input type="number" min={6} max={24} value={settings.fontSize} className={styles.num}
            onChange={e => set('fontSize', Math.max(6, Number(e.target.value) || 11))} /></label>
          {settings.layout === 'rectangular' && (
            <label>Row <input type="number" min={4} max={40} value={settings.rowHeight} className={styles.num}
              onChange={e => set('rowHeight', Math.max(4, Number(e.target.value) || 14))} /></label>
          )}
          <label>{settings.layout === 'rectangular' ? 'Width' : 'Diameter'} <input type="number" min={100} max={5000} step={50}
            value={settings.size} className={styles.num} onChange={e => set('size', Math.max(100, Number(e.target.value) || 700))} /></label>
          <span className={styles.group}>
            <button className="btn" aria-label="Zoom out" onClick={() => set('zoom', Math.max(0.1, +(settings.zoom / 1.25).toFixed(3)))}>−</button>
            <button className="btn" title="Reset zoom" onClick={() => set('zoom', 1)}>{Math.round(settings.zoom * 100)}%</button>
            <button className="btn" aria-label="Zoom in" onClick={() => set('zoom', Math.min(10, +(settings.zoom * 1.25).toFixed(3)))}>+</button>
          </span>
        </div>
        <div className={styles.row}>
          <label><input type="checkbox" checked={settings.showInternalNames} onChange={e => set('showInternalNames', e.target.checked)} /> Clade names</label>
          <label>Collapsed clades <select value={settings.wedge} onChange={e => set('wedge', e.target.value as TreeSettings['wedge'])}>
            <option value="isosceles">Isosceles</option>
            <option value="scalene">Scalene</option>
          </select></label>
          <label><input type="checkbox" checked={settings.collapsedCounts}
            onChange={e => set('collapsedCounts', e.target.checked)} /> Tip counts</label>
          {settings.layout !== 'unrooted' && <label>Tip labels <select value={settings.tipLabels === 'tips' ? 'tips' : settings.tipJustify}
            onChange={e => {
              const v = e.target.value
              setSettings(s => v === 'tips' ? { ...s, tipLabels: 'tips' } : { ...s, tipLabels: 'aligned', tipJustify: v as TreeSettings['tipJustify'] })
            }}>
            <option value="tips">at tips</option>
            <option value="left">aligned, left-justified</option>
            <option value="right">aligned, right-justified</option>
          </select></label>}
          <label><input type="checkbox" checked={settings.legend} onChange={e => set('legend', e.target.checked)} /> Size legend</label>
        </div>
        <div className={styles.row}>
          <button className="btn" onClick={() => { const op = midpointOp(ix); if (op) push(op) }}>Midpoint root</button>
          <button className="btn" disabled={!state.root} onClick={() => push({ op: 'reroot', split: null, at: 0 })}>Restore file rooting</button>
          <button className="btn" disabled={!state.collapsed.size} onClick={() => push({ op: 'expand_all' })}>Expand all</button>
          <button className="btn" onClick={() => setImportOpen(o => !o)}>Import from…</button>
          <span className={styles.group}>
            <button className="btn" disabled={!ops.length} onClick={undo} title="Ctrl+Z">Undo</button>
            <button className="btn" disabled={!redo.length} onClick={redoOne} title="Ctrl+Shift+Z">Redo</button>
          </span>
          <input className={styles.search} placeholder="Find tip" aria-label="Find tip" value={search}
            onChange={e => { setSearch(e.target.value); setMatchIdx(-1) }}
            onKeyDown={e => { if (e.key === 'Enter') jumpMatch(e.shiftKey ? -1 : 1) }} />
          {search.trim() && <>
            <span className={styles.muted}>
              {!matches.length ? 'none' : matchIdx < 0 ? `${matches.length} found` : `${matchIdx + 1} / ${matches.length}`}
            </span>
            <span className={styles.group}>
              <button className="btn" aria-label="Previous match" disabled={!matches.length} onClick={() => jumpMatch(-1)}>‹</button>
              <button className="btn" aria-label="Next match" disabled={!matches.length} onClick={() => jumpMatch(1)}>›</button>
            </span>
          </>}
          <span className={styles.spacer} />
          <span className={styles.label}>Export</span>
          <button className="btn" onClick={exportSvg}>SVG</button>
          <button className="btn" onClick={exportPng}>PNG</button>
          <button className="btn" onClick={exportNewick}>Newick</button>
          <select value={settings.newickInternal} title="Internal node labels in the Newick export"
            onChange={e => set('newickInternal', e.target.value as TreeSettings['newickInternal'])}>
            <option value="support">with support</option>
            <option value="names">with clade names</option>
            <option value="both">with names and support (NHX)</option>
          </select>
          {source.study && <AddToReport study={source.study} kind="tree" className="btn" defaultTitle={baseName}
            make={async () => ({ blob: new Blob([svgMarkup()], { type: 'image/svg+xml' }), ext: '.svg' })} />}
        </div>
      </div>

      <div className={styles.body}>
        <div className={styles.canvas} ref={scrollRef} onScroll={onScroll}>
          <svg ref={svgRef} width={layout.width} height={layout.height + legendH} className={styles.svg}
            fontFamily="Arial, Helvetica, sans-serif" fontSize={settings.fontSize}
            onClick={e => { if (e.target === svgRef.current) setSelClade(null) }}>
            {legendRows.length > 0 && <Legend rows={legendRows} settings={settings} />}
            <g transform={legendH ? `translate(0 ${legendH})` : undefined}>
              <TreeShapes layout={layout} settings={settings} selSubtree={selSubtree} matchSet={matchSet}
                onSelect={d => setSelClade(d.clade)}
                onToggle={toggleCollapse}
                onRename={startRename} />
              {layout.symbols.map((sy, k) => {
                const st = settings.symbols[sy.dataset]
                const title = sy.dataset === 'support'
                  ? `Support ${sy.node.support}`
                  : `${fmt(sy.value, 3)} placement${sy.value === 1 ? '' : 's'}\n` + sy.queries
                      .map(q => `${src.pqueries[q.query].names.join(', ')} (LWR ${fmt(q.lwr, 3)})`).join('\n')
                return (
                  <path key={k} d={symbolPath(st.shape, sy.x, sy.y, sy.size)}
                    fill={st.fill.hex} fillOpacity={st.fill.alpha}
                    stroke={st.stroke.hex} strokeOpacity={st.stroke.alpha} strokeWidth={1}
                    className={styles.placement}
                    onClick={e => {
                      e.stopPropagation()
                      setSelClade(sy.node.clade)
                      if (sy.queries.length === 1) setSelQuery(sy.queries[0].query)
                    }}>
                    <title>{title}</title>
                  </path>
                )
              })}
              {selQuery !== null && layout.placements.filter(p => p.query === selQuery).map((p, k) => (
                <circle key={`sq${k}`} data-ui="" cx={p.x} cy={p.y} r={4 + 4 * Math.sqrt(p.lwr)}
                  fill="none" stroke={PLACE_SELECTED} strokeWidth={2} pointerEvents="none" />
              ))}
              {settings.placementLabels && layout.placements.filter(p => p.rank === 0).map((p, k) => (
                <text key={`pl${k}`} x={p.x + 5} y={p.y - 6} fontSize={settings.fontSize * 0.85}
                  fill={settings.symbols.placements.fill.hex}>
                  {src.pqueries[p.query].names.join(', ')}
                </text>
              ))}
            </g>
          </svg>
          {editing && (
            <input autoFocus className={styles.renameInput}
              style={{ left: editing.left, top: editing.top, width: editing.width, fontSize: settings.fontSize + 1 }}
              value={editing.value}
              onChange={e => setEditing({ ...editing, value: e.target.value })}
              onKeyDown={e => {
                if (e.key === 'Enter') commitRename()
                else if (e.key === 'Escape') setEditing(null)
              }}
              onBlur={commitRename} />
          )}
        </div>

        <aside className={styles.side}>
          {selected ? (
            <Inspector
              d={selected} ix={ix} names={state.names}
              placements={selectedPlacements.map(p => ({ name: src.pqueries[p.query].names.join(', '), lwr: p.lwr, rank: p.rank, query: p.query }))}
              onRename={(key, name) => push(renameOp(key, name, selected.clade))}
              onStyle={style => push({ op: 'style', clade: selected.clade, style })}
              onRotate={() => push({ op: 'rotate', clade: selected.clade })}
              onToggle={() => toggleCollapse(selected)}
              onReroot={() => reroot(selected)}
              onParent={() => selected.parent && setSelClade(selected.parent.clade)}
              onQuery={selectQuery}
            />
          ) : (
            <div className="card">
              <div className="card-title">Selection</div>
              <p className={styles.muted}>
                Click a branch to select it. Keys: C collapse, R reroot, O rotate, F2 rename, ↑ parent, Ctrl+Z undo.
              </p>
            </div>
          )}
          <div className="card">
            <div className="card-title">Annotations</div>
            {src.pqueries.length > 0 && (
              <DatasetControls title="Placements" style={settings.symbols.placements}
                onChange={p => setSym('placements', p)}>
                <label>Count <select value={settings.placements}
                  onChange={e => set('placements', e.target.value as TreeSettings['placements'])}>
                  <option value="best">best placement</option>
                  <option value="all">all candidates, LWR-weighted</option>
                </select></label>
                {settings.placements === 'all' && (
                  <label>LWR at least <input type="number" min={0} max={1} step={0.05} value={settings.placementMinLwr} className={styles.num}
                    onChange={e => set('placementMinLwr', Math.min(1, Math.max(0, Number(e.target.value) || 0)))} /></label>
                )}
                <label><input type="checkbox" checked={settings.placementLabels}
                  onChange={e => set('placementLabels', e.target.checked)} /> Label queries</label>
              </DatasetControls>
            )}
            {support.size > 0 ? (
              <DatasetControls title="Support" style={settings.symbols.support}
                onChange={p => setSym('support', p)}>
                {state.supportFrom && (
                  <span className={styles.muted}>
                    {state.support.size} values from {state.supportFrom}{' '}
                    <button className={styles.linkBtn} onClick={() => push({ op: 'support', from: null, values: [] })}>Remove</button>
                  </span>
                )}
                <label><input type="checkbox" checked={settings.showSupport}
                  onChange={e => set('showSupport', e.target.checked)} /> Values as text</label>
                <label>At least <input type="number" min={0} step="any" value={settings.supportMin} className={styles.num}
                  onChange={e => set('supportMin', Number(e.target.value) || 0)} /></label>
              </DatasetControls>
            ) : (
              <p className={styles.muted}>No support values. Use Import from… to take them from another tree.</p>
            )}
          </div>
          {src.pqueries.length > 0 && (
            <div className="card">
              <div className="card-title">Placements</div>
              <div className={styles.queryList}>
                <table>
                  <thead><tr><th>Query</th><th>LWR</th><th>Location</th></tr></thead>
                  <tbody>
                    {queryRows.map(r => (
                      <tr key={r.qi} className={selQuery === r.qi ? styles.rowOn : undefined} onClick={() => selectQuery(r.qi)}
                        tabIndex={0} aria-selected={selQuery === r.qi}
                        onKeyDown={e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); selectQuery(r.qi) } }}>
                        <td>{r.name}</td>
                        <td title={`${r.n} candidate${r.n === 1 ? '' : 's'}`}>{fmt(r.lwr, 3)}{r.n > 1 ? ` (${r.n})` : ''}</td>
                        <td className={styles.where}>{r.where}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </aside>
      </div>
    </>
  )
}
