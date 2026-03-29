import { useEffect, useRef } from 'react'
import Plotly, { type Data, type Layout } from 'plotly.js-dist-min'

interface PlotlySpec {
  data:   Data[]
  layout: Partial<Layout>
}

interface Props {
  figure: unknown
  className?: string
}

export function PlotlyChart({ figure, className }: Props) {
  const ref = useRef<HTMLDivElement>(null)

  useEffect(() => {
    if (!ref.current || !figure) return
    const spec = figure as PlotlySpec
    Plotly.react(ref.current, spec.data ?? [], {
      autosize: true,
      height: ref.current.clientWidth * 0.6,
      margin: { l: 60, r: 30, t: 40, b: 50 },
      ...spec.layout,
    }, { responsive: true, displaylogo: false })
    return () => { ref.current && Plotly.purge(ref.current) }
  }, [figure])

  return <div ref={ref} className={className} style={{ width: '100%', aspectRatio: '5 / 3' }} />
}
