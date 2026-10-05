// (c) 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useRef, useState } from 'react'
import { api } from '../api/client'
import { useAnalysisSource } from './analysisSource'
import { errorMessage } from '../api/errorMessage'
import { AnalysisChart } from './AnalysisChart'
import { useToast } from './Toast'
import { FACET_DIMENSIONS, type CategorySet, type ComparisonRunSpec, type FacetDimension } from '../api/types'

//## Select / checkbox style tokens (matching DiversityPanel.tsx)
const SELECT_STYLE: React.CSSProperties = {
  padding: '4px 8px',
  borderRadius: 4,
  border: '1px solid var(--color-border)',
  fontSize: '.82rem',
  background: 'var(--color-bg)',
}

const LABEL_STYLE: React.CSSProperties = {
  fontSize: '.82rem',
  display: 'flex',
  alignItems: 'center',
  gap: 4,
}

export function TaxaCompositionChart({
  study,
  run,
  group,
  runs,
  subgroups,
  defaultTag,
  table,
  subgroup: controlledSubgroup,
}: {
  study: string
  run?: string
  group?: string | null
  runs?: ComparisonRunSpec[]
  subgroups?: string[]
  defaultTag: 'rank' | 'category'
  table?: string
  // When provided (even when null), this value is used as the sub-group for
  // requests and the chart's internal sub-group selector is hidden.
  subgroup?: string | null
}) {
  const source = useAnalysisSource()
  const toast = useToast()

  //## Tag selector: 'rank' or 'category'
  const [tag, setTag] = useState<'rank' | 'category'>(defaultTag)

  //## Rank state (used when tag='rank')
  const [ranks, setRanks] = useState<string[]>([])
  const [rank, setRank] = useState<string | null>(null)
  const [topN, setTopN] = useState(15)

  //## Category-set state (used when tag='category')
  const [catSets, setCatSets] = useState<CategorySet[]>([])
  const [catSet, setCatSet] = useState<string>('default')

  //## Shared controls
  const [relative, setRelative] = useState(true)
  const [mode, setMode] = useState<'stacked' | 'grouped'>('stacked')
  const [keepEmpty, setKeepEmpty] = useState(false)
  // Grid layout: one panel per (row dimension value, column dimension value).
  const [facet, setFacet] = useState(false)
  const [facetRows, setFacetRows] = useState<FacetDimension>('run')
  const [facetCols, setFacetCols] = useState<FacetDimension>('subgroup')
  // null = All, "__pool__" = Pool, any other string = a specific sub-group.
  // When controlledSubgroup is provided, the internal selector is suppressed
  // and this state is ignored in favour of the prop.
  const [internalSubgroup, setInternalSubgroup] = useState<string | null>(null)
  const isControlled = controlledSubgroup !== undefined
  const subgroup = isControlled ? controlledSubgroup : internalSubgroup

  //## Chart state
  const [figure, setFigure] = useState<unknown>(null)
  const [loading, setLoading] = useState(false)
  // Only the response to the latest request, for the current inputs, is shown.
  const request = useRef(0)

  //## Derived flags
  const isCrossRun = runs !== undefined && runs.length > 0
  // Faceting splits the selected runs across two dimensions, so it is only on
  // offer where several runs are in play.
  const canFacet = isCrossRun
  const isFaceted = canFacet && facet
  // Show the internal sub-group selector only when not controlled externally, and
  // when subgroups has >= 2 entries or we are in cross-run mode. A facet grid
  // decides the scoping itself, so the selector would only contradict it.
  const showSubgroupSelector =
    !isControlled && !isFaceted &&
    (isCrossRun || (subgroups !== undefined && subgroups.length >= 2))
  const effectiveTable = table ?? 'merged'
  // Stable scalars for the reference run and group, shared by both single-run
  // and cross-run paths. These are primitive values so the rank-fetch effect
  // depends on them without churn from freshly-constructed runs arrays.
  const refRun = run ?? runs?.[0]?.run
  const refGroup = group ?? runs?.[0]?.group

  //## Fetch taxonomy ranks when tag='rank' or the reference run/group changes
  useEffect(() => {
    if (tag !== 'rank') return
    if (!refRun) { setRanks([]); return }
    let cancelled = false
    api.analysis
      .ranks(study, refRun, { table: effectiveTable, group: refGroup ?? undefined })
      .then(result => {
        if (cancelled) return
        setRanks(result)
        setRank(current =>
          current && result.includes(current)
            ? current
            : result[result.length - 1] ?? null
        )
      })
      .catch(() => { if (!cancelled) setRanks([]) })
    return () => { cancelled = true }
  }, [study, refRun, refGroup, tag, effectiveTable])

  //## Fetch category sets when tag='category'
  useEffect(() => {
    if (tag !== 'category') return
    api.composition.categorySets().then(sets => {
      setCatSets(sets)
      setCatSet(current =>
        sets.some(cs => cs.name === current) ? current : sets[0]?.name ?? 'default'
      )
    }).catch(() => {})
  }, [tag])

  // Panel rows in the rendered figure, read back from the layout Plotly grid.
  // A single bar chart carries no grid, hence the fallback of one row.
  const gridRows =
    (figure as { layout?: { grid?: { rows?: number } } } | null)?.layout?.grid?.rows ?? 1

  // Everything a request is built from. A figure answers the inputs it was
  // computed from, so new inputs discard it, along with any request still in
  // flight for the old ones. Keyed on content so a new-but-equal `runs` array
  // does not clear the chart.
  const inputsKey = JSON.stringify([
    study, run, group, runs, effectiveTable, source, tag, rank, catSet, topN,
    relative, mode, keepEmpty, isFaceted, facetRows, facetCols, subgroup,
  ])
  useEffect(() => {
    request.current++
    setFigure(null)
    setLoading(false)
  }, [inputsKey])

  //## Compute handler
  /** Build the chart for the current inputs and show it, unless the inputs changed meanwhile. */
  const handleShow = async () => {
    const value = tag === 'rank' ? rank : catSet
    if (!value) return
    // Guard: a single-run chart requires a run identifier.
    if (!isCrossRun && run === undefined) {
      toast.error('No run specified')
      return
    }
    const req = ++request.current
    setLoading(true)
    try {
      const body = {
        table: effectiveTable,
        tag,
        value,
        relative,
        mode,
        keep_empty: keepEmpty,
        ...(tag === 'rank' ? { top_n: topN } : {}),
      }
      let fig: unknown
      if (isFaceted) {
        fig = await api.analysis.chartFacet(study, {
          ...body,
          runs: runs!,
          rows: facetRows,
          cols: facetCols,
        })
      } else if (isCrossRun) {
        fig = await api.analysis.chartCompare(study, {
          ...body,
          subgroup: subgroup ?? null,
          runs: runs!,
        })
      } else {
        fig = await api.analysis.chart(study, run!,
          { ...body, subgroup: subgroup ?? null, ...(source ? { source } : {}) }, group)
      }
      if (req === request.current) setFigure(fig)
    } catch (err) {
      if (req !== request.current) return
      setFigure(null)
      toast.error(`Chart failed: ${errorMessage(err)}`)
    } finally {
      if (req === request.current) setLoading(false)
    }
  }

  return (
    <div>
      {/* Control strip */}
      <div style={{ display: 'flex', gap: 8, alignItems: 'center', flexWrap: 'wrap', marginBottom: 12 }}>
        {/* Tag-by selector */}
        <label style={LABEL_STYLE}>
          Tag by:
          <select
            value={tag}
            onChange={e => setTag(e.target.value as 'rank' | 'category')}
            style={SELECT_STYLE}
          >
            <option value="rank">Rank</option>
            <option value="category">Category</option>
          </select>
        </label>

        {/* Value selector: ranks when tag='rank', category sets when tag='category' */}
        {tag === 'rank' ? (
          <select
            aria-label="Rank"
            value={rank ?? ''}
            onChange={e => setRank(e.target.value)}
            style={SELECT_STYLE}
          >
            {ranks.map(r => <option key={r} value={r}>{r}</option>)}
          </select>
        ) : (
          <select
            aria-label="Category set"
            value={catSet}
            onChange={e => setCatSet(e.target.value)}
            style={SELECT_STYLE}
          >
            {catSets.map(cs => <option key={cs.name} value={cs.name}>{cs.label}</option>)}
          </select>
        )}

        {/* top_n: shown only when tag='rank' */}
        {tag === 'rank' && (
          <label style={LABEL_STYLE}>
            Top N:
            <input
              type="number"
              min={1}
              max={100}
              value={topN}
              onChange={e => {
                const n = parseInt(e.target.value, 10)
                setTopN(Number.isNaN(n) ? 15 : n)
              }}
              style={{
                width: 52,
                padding: '3px 6px',
                borderRadius: 4,
                border: '1px solid var(--color-border)',
                fontSize: '.82rem',
              }}
            />
          </label>
        )}

        {/* Keep zero-read samples as blank slots on the axis */}
        <label style={LABEL_STYLE} title="Show samples with no reads as blank slots on the axis">
          <input
            type="checkbox"
            checked={keepEmpty}
            onChange={e => setKeepEmpty(e.target.checked)}
          />
          Show empty
        </label>

        {/* Relative checkbox */}
        <label style={LABEL_STYLE}>
          <input
            type="checkbox"
            checked={relative}
            onChange={e => setRelative(e.target.checked)}
          />
          Relative
        </label>

        {/* Stacked / Grouped selector */}
        <select
          aria-label="Bar mode"
          value={mode}
          onChange={e => setMode(e.target.value as 'stacked' | 'grouped')}
          style={SELECT_STYLE}
        >
          <option value="stacked">Stacked</option>
          <option value="grouped">Grouped</option>
        </select>

        {/* Facet grid: one panel per (row value, column value) pair */}
        {canFacet && (
          <label style={LABEL_STYLE} title="Split the selected runs into a grid of panels">
            <input
              type="checkbox"
              checked={facet}
              onChange={e => setFacet(e.target.checked)}
            />
            Grid
          </label>
        )}
        {isFaceted && (
          <>
            <label style={LABEL_STYLE}>
              Rows:
              <select
                value={facetRows}
                onChange={e => {
                  const next = e.target.value as FacetDimension
                  setFacetRows(next)
                  // The two axes must name different dimensions, so bump the other one.
                  if (next === facetCols) {
                    setFacetCols(FACET_DIMENSIONS.find(d => d !== next)!)
                  }
                }}
                style={SELECT_STYLE}
              >
                {FACET_DIMENSIONS.map(d => <option key={d} value={d}>{d}</option>)}
              </select>
            </label>
            <label style={LABEL_STYLE}>
              Columns:
              <select
                value={facetCols}
                onChange={e => {
                  const next = e.target.value as FacetDimension
                  setFacetCols(next)
                  if (next === facetRows) {
                    setFacetRows(FACET_DIMENSIONS.find(d => d !== next)!)
                  }
                }}
                style={SELECT_STYLE}
              >
                {FACET_DIMENSIONS.map(d => <option key={d} value={d}>{d}</option>)}
              </select>
            </label>
          </>
        )}

        {/* Internal sub-group selector: hidden when the prop controls the value */}
        {showSubgroupSelector && (
          <select
            value={internalSubgroup ?? ''}
            onChange={e => setInternalSubgroup(e.target.value || null)}
            style={SELECT_STYLE}
          >
            <option value="">All</option>
            {(subgroups ?? []).map(sg => (
              <option key={sg} value={sg}>{sg}</option>
            ))}
            <option value="__pool__">Pool</option>
          </select>
        )}

        {/* Compute button */}
        <button
          className="btn btn-primary"
          onClick={handleShow}
          disabled={loading || (tag === 'rank' ? !rank : !catSet)}
        >
          {loading ? 'Computing…' : 'Show'}
        </button>
      </div>

      {/* A grid needs headroom for its rows. */}
      {figure != null && (
        <AnalysisChart
          study={study}
          figure={figure}
          heightRatio={Math.min(0.3 * gridRows, 0.9)}
        />
      )}
    </div>
  )
}
