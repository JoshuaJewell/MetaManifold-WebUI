// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import Plotly from 'plotly.js-dist-min'
import { AddToReport } from './AddToReport'
import { PlotlyChart } from './PlotlyChart'

const figureTitle = (fig: unknown): string | null => {
  const t = (fig as { layout?: { title?: unknown } })?.layout?.title
  const text = typeof t === 'string' ? t : (t as { text?: string } | undefined)?.text
  return text ? text.replace(/<[^>]+>/g, '') : null
}

// The chart as standalone SVG, at a fixed print width.
async function figureSvg(fig: unknown, heightRatio: number): Promise<{ blob: Blob; ext: string }> {
  const f = fig as { data?: unknown[]; layout?: Record<string, unknown> }
  const width = 1200
  const url = await Plotly.toImage(
    { data: (f.data ?? []) as Record<string, unknown>[], layout: { ...(f.layout ?? {}), width, height: Math.round(width * heightRatio) } },
    { format: 'svg', width, height: Math.round(width * heightRatio) })
  const svg = decodeURIComponent(url.slice(url.indexOf(',') + 1))
  return { blob: new Blob([svg], { type: 'image/svg+xml' }), ext: '.svg' }
}

/** An analysis chart with its Add to report button. */
export function AnalysisChart({ study, figure, heightRatio }: {
  study: string
  figure: unknown
  heightRatio?: number
}) {
  return (
    <div>
      <div style={{ display: 'flex', justifyContent: 'flex-end', gap: 6, marginBottom: 4 }}>
        <AddToReport study={study} kind="figure" className="btn" defaultTitle={figureTitle(figure) ?? 'Chart'}
          make={() => figureSvg(figure, heightRatio ?? 0.6)} />
      </div>
      <PlotlyChart figure={figure} heightRatio={heightRatio} />
    </div>
  )
}
