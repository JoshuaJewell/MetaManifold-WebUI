// SPDX-License-Identifier: AGPL-3.0-only
import { useEffect, useRef, useState } from 'react'
import type { Data, Layout } from 'plotly.js-dist-min'

interface PlotlySpec {
  data:   Data[]
  layout: Partial<Layout>
}

interface Props {
  figure: unknown
  className?: string | undefined
  /** Height as a ratio of container width (default 0.6). */
  heightRatio?: number | undefined
}

// Lazy-load Plotly to keep it out of main bundle — isolated to plotly chunk via vite.config.ts manualChunks
// Future gossamer UI (ui/) will render charts server-side or via lightweight client without 4.6MB bundle
export function PlotlyChart({ figure, className, heightRatio = 0.6 }: Props) {
  const wrapRef = useRef<HTMLDivElement>(null)
  const plotRef = useRef<HTMLDivElement>(null)
  const ready   = useRef(false)
  const hasInitialized = useRef(false)
  const [dims, setDims] = useState<{ w: number; h: number } | null>(null)
  const [Plotly, setPlotly] = useState<any>(null)
  const [loadError, setLoadError] = useState<string | null>(null)

  // Dynamic import of plotly — only loaded when chart is rendered
  useEffect(() => {
    let cancelled = false
    import('plotly.js-dist-min')
      .then(mod => {
        if (!cancelled) setPlotly(mod.default)
      })
      .catch(e => {
        if (!cancelled) setLoadError(e instanceof Error ? e.message : String(e))
      })
    return () => { cancelled = true }
  }, [])

  useEffect(() => {
    const wrap = wrapRef.current
    const plot = plotRef.current
    if (!wrap || !plot || !figure || !Plotly) return
    const spec = figure as PlotlySpec
    const data = (spec.data ?? []).map((t: Record<string, unknown>) =>
      t['type'] === 'box' && t['width'] == null ? { ...t, width: 0.9 } : t
    )
    Plotly.react(plot, data as Data[], {
      autosize: true,
      height: wrap.clientHeight,
      width:  wrap.clientWidth,
      margin: { l: 60, r: 30, t: 40, b: 50 },
      ...spec.layout,
      font:   { size: 16, ...(spec.layout?.['font']   as object | undefined) },
      legend: { font: { size: 20 }, ...(spec.layout?.['legend'] as object | undefined) },
    }, { responsive: true, displaylogo: false })
    ready.current = true

    if (!hasInitialized.current && wrap.clientWidth > 0 && wrap.clientHeight > 0) {
      setDims({ w: wrap.clientWidth, h: wrap.clientHeight })
      hasInitialized.current = true
    }

    return () => { ready.current = false; plot && Plotly.purge(plot) }
  }, [figure, heightRatio, Plotly])

  useEffect(() => {
    const wrap = wrapRef.current
    const plot = plotRef.current
    if (!wrap || !plot || !Plotly) return
    const ro = new ResizeObserver(() => {
      if (!ready.current) return
      const w = wrap.clientWidth
      const h = wrap.clientHeight
      Plotly.relayout(plot, { height: h, width: w })
      setDims(prev => (prev?.w === w && prev?.h === h) ? prev : { w, h })
    })
    ro.observe(wrap)
    return () => ro.disconnect()
  }, [Plotly])

  function applyDims(w: number, h: number) {
    if (!Number.isFinite(w) || !Number.isFinite(h) || w < 1 || h < 1) return
    const wrap = wrapRef.current
    if (!wrap) return
    wrap.style.width  = `${w}px`
    wrap.style.height = `${h}px`
    setDims({ w, h })
  }

  const inputStyle: React.CSSProperties = {
    width: 70, marginLeft: 4, padding: '2px 4px',
    borderRadius: 3, border: '1px solid var(--color-border)',
    fontSize: 'inherit', background: 'var(--color-bg)',
    color: 'var(--color-fg)',
  }

  if (loadError) {
    return <div style={{ padding: 12, color: 'var(--color-danger, #c00)' }}>Plotly failed to load: {loadError}</div>
  }

  if (!Plotly) {
    return <div style={{ padding: 12 }}>Loading chart...</div>
  }

  return (
    <div style={{ overflowX: 'auto' }}>
      <div
        ref={wrapRef}
        className={className}
        style={{
          width: '100%',
          aspectRatio: `${1 / heightRatio}`,
          resize: 'both',
          overflow: 'hidden',
          minHeight: 120,
          minWidth: 300,
        }}
      >
        <div ref={plotRef} style={{ width: '100%', height: '100%' }} />
      </div>
      {dims && (
        <div style={{ display: 'flex', gap: 10, alignItems: 'center', marginTop: 4,
                      fontSize: '.75rem', color: 'var(--color-muted-fg)' }}>
          <label>
            W (px)
            <input
              type="number" min={300} step={10} value={dims.w}
              onChange={e => applyDims(Number(e.target.value), dims.h)}
              style={inputStyle}
            />
          </label>
          <label>
            H (px)
            <input
              type="number" min={120} step={10} value={dims.h}
              onChange={e => applyDims(dims.w, Number(e.target.value))}
              style={inputStyle}
            />
          </label>
        </div>
      )}
    </div>
  )
}
