// SPDX-License-Identifier: AGPL-3.0-only
// (c) 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).

// ChartEditorInner — isolated heavy editor with eval
// This file imports react-chart-editor which uses eval internally (100+ Vite warnings).
// It is intentionally isolated via React.lazy() in ChartCustomiser.tsx and via manualChunks
// in vite.config.ts (chunk 'chartEditor'). It is only loaded when user clicks "Customise".
// Eval is sandboxed to this chunk, not in main bundle.
// Future: gossamer (Stipple/Vue UI in ui/) will replace this entirely — no JS build, no eval,
// server-rendered via Genie/Stipple/StippleUI. See docs/migration/gossamer-vite-hybrid.md.

import { useEffect, useState } from 'react'
import type { ChartEditorState, ChartEditorUpdateHandler } from '../types/components'

// CSS is imported here but will be code-split into chartEditor chunk
import 'react-chart-editor/lib/react-chart-editor.css'
import './chartEditor.css'

// Dynamic imports to avoid top-level eval execution
// Plotly and react-chart-editor are both heavy and contain eval
type PlotlyType = typeof import('plotly.js-dist-min').default
type EditorModule = typeof import('react-chart-editor')

export default function ChartEditorInner({ state, onUpdate }: {
  state: ChartEditorState
  onUpdate: ChartEditorUpdateHandler
}) {
  const [Plotly, setPlotly] = useState<PlotlyType | null>(null)
  const [Editor, setEditor] = useState<EditorModule | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    let cancelled = false
    async function load() {
      try {
        // Check feature flag — allows strict CSP to disable editor
        // @ts-ignore defined in vite.config.ts
        if (typeof __ENABLE_CHART_EDITOR__ !== 'undefined' && !__ENABLE_CHART_EDITOR__) {
          throw new Error('Chart editor disabled by feature flag __ENABLE_CHART_EDITOR__ — use gossamer UI at /studies')
        }
        const [plotlyMod, editorMod] = await Promise.all([
          import('plotly.js-dist-min'),
          import('react-chart-editor')
        ])
        if (!cancelled) {
          setPlotly(plotlyMod.default as unknown as PlotlyType)
          setEditor(editorMod as unknown as EditorModule)
        }
      } catch (e) {
        if (!cancelled) {
          setError(e instanceof Error ? e.message : String(e))
        }
      }
    }
    load()
    return () => { cancelled = true }
  }, [])

  if (error) {
    return (
      <div style={{ padding: 24, border: '1px solid var(--color-danger, #c00)', borderRadius: 6 }}>
        <h3 style={{ margin: '0 0 8px' }}>Chart editor unavailable</h3>
        <p style={{ margin: '0 0 8px', fontSize: '.9rem' }}>{error}</p>
        <p style={{ margin: 0, fontSize: '.8rem', color: 'var(--color-muted-fg)' }}>
          This editor uses eval and is isolated to its own chunk. In strict CSP environments,
          disable it via <code>__ENABLE_CHART_EDITOR__=false</code> in vite.config.ts and use
          the gossamer UI at <code>/studies</code> (Stipple/Vue, no JS build, no eval).
          See docs/migration/gossamer-vite-hybrid.md
        </p>
      </div>
    )
  }

  if (!Plotly || !Editor) {
    return <div style={{ padding: 24 }}>Loading editor (plotly + chart editor)...</div>
  }

  const PlotlyEditor = Editor.default
  const { PanelMenuWrapper, StyleLayoutPanel, StyleTracesPanel, StyleAxesPanel, StyleLegendPanel } = Editor

  return (
    <PlotlyEditor
      data={state.data} layout={state.layout} frames={state.frames}
      config={{ editable: true, displaylogo: false }} plotly={Plotly}
      onUpdate={onUpdate} useResizeHandler
    >
      <PanelMenuWrapper>
        <StyleLayoutPanel group="Style" name="General" />
        <StyleTracesPanel group="Style" name="Traces" />
        <StyleAxesPanel group="Style" name="Axes" />
        <StyleLegendPanel group="Style" name="Legend" />
      </PanelMenuWrapper>
    </PlotlyEditor>
  )
}
