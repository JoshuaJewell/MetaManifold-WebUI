// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useCallback, useEffect, useMemo, useState } from 'react'
import { api } from '../api/client'
import { useAnalysisSource } from './analysisSource'
import { errorMessage } from '../api/errorMessage'
import { DataTable } from './DataTable'
import { NameDialog } from './NameDialog'
import { useToast } from './Toast'
import { discoverResultsTables, type AnalysisOption } from './annotationShared'
import type {
  AnnotationSource, CategorySet, ColFilter, CompositionBuildResult,
  TableQuery,
} from '../api/types'

// Fixed colour of the catch-all "Unassigned" bucket, mirroring the backend.
const UNASSIGNED_COLOUR = '#95a5a6'

export function CompositionPanel({
  study, run, group, subgroups, source: configuredSource,
}: {
  study: string
  run: string
  group?: string
  subgroups?: string[]
  // Pipeline-configured annotation source, resolved from tagging.source in the
  // run config cascade. Passed by RunView; defaults to VSEARCH when absent.
  source?: AnnotationSource
}) {
  // The analysis tab's classifier wins over the run's configured one.
  const source = useAnalysisSource() ?? configuredSource
  const toast = useToast()
  const effectiveSource: AnnotationSource = source ?? 'VSEARCH'

  //## Category set selection
  const [categorySets, setCategorySets] = useState<CategorySet[]>([])
  const [selectedCatSet, setSelectedCatSet] = useState<string>('default')

  //## Sub-group scope: null = All, or a specific sub-group name
  const [subgroup, setSubgroup] = useState<string | null>(null)

  //## Results table and breakdown, with the same choices as the composition chart
  const [tableOptions, setTableOptions] = useState<AnalysisOption[]>([])
  const [table, setTable] = useState('merged')
  const [tag, setTag] = useState<'rank' | 'category'>('category')
  const [ranks, setRanks] = useState<string[]>([])
  const [rank, setRank] = useState<string | null>(null)

  useEffect(() => {
    let cancelled = false
    discoverResultsTables(study, run, group ?? null).then(opts => {
      if (cancelled) return
      setTableOptions(opts)
      setTable(current => opts.some(o => o.key === current) ? current
        : opts.some(o => o.key === 'merged') ? 'merged' : opts[0]?.key ?? current)
    })
    return () => { cancelled = true }
  }, [study, run, group])

  useEffect(() => {
    if (tag !== 'rank') return
    let cancelled = false
    api.analysis.ranks(study, run, { table, group: group ?? null, ...(source ? { source } : {}) })
      .then(rs => {
        if (cancelled) return
        setRanks(rs)
        // Deepest rank by default, matching the chart route's default.
        setRank(current => current && rs.includes(current) ? current : rs[rs.length - 1] ?? null)
      })
      .catch(() => { if (!cancelled) setRanks([]) })
    return () => { cancelled = true }
  }, [study, run, group, table, tag])

  const byCategory = tag === 'category'

  //## Summary state
  const [summary, setSummary] = useState<CompositionBuildResult | null>(null)

  //## Live colour editing
  // Pending per-category colour overrides, keyed by category name. Reset
  // whenever the selected set changes so edits never leak across sets.
  const [colourOverrides, setColourOverrides] = useState<Record<string, string>>({})
  // Guards the Save / Save as / Delete buttons against concurrent requests.
  const [catSetBusy, setCatSetBusy] = useState(false)
  const [showSaveAs, setShowSaveAs] = useState(false)
  useEffect(() => {
    setColourOverrides({})
  }, [selectedCatSet])

  const activeCatSet = useMemo(
    () => categorySets.find(cs => cs.name === selectedCatSet) ?? null,
    [categorySets, selectedCatSet],
  )

  const categoryColourMap = useMemo(() => {
    const map: Record<string, string> = {
      Unassigned: activeCatSet?.unassigned_colour ?? UNASSIGNED_COLOUR,
    }
    for (const c of activeCatSet?.categories ?? []) {
      if (c.colour) map[c.name] = c.colour
    }
    return map
  }, [activeCatSet])

  // Effective colour map for rendering: the set's colours overlaid with any
  // pending overrides, the override winning.
  const effectiveColourMap = useMemo(
    () => ({ ...categoryColourMap, ...colourOverrides }),
    [categoryColourMap, colourOverrides],
  )

  // Returning a category to its set colour drops the override, so Save
  // reflects only real changes.
  const setCategoryColour = useCallback((name: string, colour: string) => {
    setColourOverrides(prev => {
      const base = categoryColourMap[name]
      const next = { ...prev }
      if (base && base.toLowerCase() === colour.toLowerCase()) {
        delete next[name]
      } else {
        next[name] = colour
      }
      return next
    })
  }, [categoryColourMap])

  const hasPendingColours = Object.keys(colourOverrides).length > 0

  // Categories that can be persisted: those declared in the selected set, plus the
  // Unassigned catch-all, whose colour the backend routes to unassigned_colour.
  const savableCategories = useMemo(
    () => new Set<string>([
      ...(activeCatSet?.categories ?? []).map(c => c.name),
      'Unassigned',
    ]),
    [activeCatSet],
  )

  // Refetch the category sets after a mutation; returns the fresh list.
  const refetchCatSets = useCallback(async () => {
    const sets = await api.composition.categorySets()
    setCategorySets(sets)
    return sets
  }, [])

  //## Category-set mutations
  const handleSaveCatSet = useCallback(async () => {
    setCatSetBusy(true)
    try {
      await api.composition.saveCategorySet(selectedCatSet, {
        base: selectedCatSet, colours: colourOverrides,
      })
      await refetchCatSets()
      setColourOverrides({})
      toast.success('Colours saved')
    } catch (err) {
      toast.error(errorMessage(err, 'Failed to save colours'))
    } finally {
      setCatSetBusy(false)
    }
  }, [selectedCatSet, colourOverrides, refetchCatSets, toast])

  const handleSaveAsCatSet = useCallback(async (newName: string) => {
    setCatSetBusy(true)
    try {
      await api.composition.saveCategorySet(newName, {
        base: selectedCatSet, colours: colourOverrides,
      })
      await refetchCatSets()
      setSelectedCatSet(newName)
      setShowSaveAs(false)
      toast.success(`Saved as ${newName}`)
    } catch (err) {
      // Rethrow so NameDialog surfaces the message inline.
      throw err
    } finally {
      setCatSetBusy(false)
    }
  }, [selectedCatSet, colourOverrides, refetchCatSets])

  const handleDeleteCatSet = useCallback(async () => {
    if (!window.confirm(`Delete category set "${selectedCatSet}"?`)) return
    setCatSetBusy(true)
    try {
      await api.composition.deleteCategorySet(selectedCatSet)
      await refetchCatSets()
      setSelectedCatSet('default')
      toast.success('Category set deleted')
    } catch (err) {
      toast.error(errorMessage(err, 'Failed to delete category set'))
    } finally {
      setCatSetBusy(false)
    }
  }, [selectedCatSet, refetchCatSets, toast])

  // Fetch category sets on mount
  useEffect(() => {
    api.composition.categorySets().then(setCategorySets).catch(() => {})
  }, [])

  //## Fetch the summary whenever its scope or breakdown changes
  useEffect(() => {
    setSummary(null)
    if (!byCategory && !rank) return
    let cancelled = false
    api.composition.summary(study, run, {
      category_set: selectedCatSet, subgroup, table, tag,
      value: byCategory ? selectedCatSet : rank!,
      ...(source ? { source } : {}),
    }, group ?? null)
      .then(r => { if (!cancelled) setSummary(r) })
      .catch(() => { if (!cancelled) setSummary(null) })
    return () => { cancelled = true }
  }, [study, run, group, selectedCatSet, subgroup, table, tag, rank, byCategory, source])

  //## Paginated table fetchers (include selectedCatSet so the backend tags the right column)
  const fetcher = useCallback(
    (q: TableQuery) =>
      api.composition.query(study, run, effectiveSource, q, group, selectedCatSet, table),
    [study, run, effectiveSource, group, selectedCatSet, table],
  )

  const distinctFetcher = useCallback(
    (column: string, activeFilters?: Record<string, ColFilter>, keywordFilter?: string) =>
      api.composition.distinct(study, run, effectiveSource, column, activeFilters,
        group, selectedCatSet, keywordFilter, table),
    [study, run, effectiveSource, group, selectedCatSet, table],
  )

  return (
    <div className="card">
      {/* Table, breakdown and sub-group selectors */}
      <div style={{ display: 'flex', gap: 8, alignItems: 'center', marginBottom: 12, flexWrap: 'wrap' }}>
        <label style={{ fontWeight: 600, fontSize: '.85rem' }}>Table:</label>
        <select
          value={table}
          onChange={e => setTable(e.target.value)}
          style={{ padding: '4px 8px', borderRadius: 4, border: '1px solid var(--color-border)' }}
        >
          {tableOptions.length === 0 && <option value={table}>{table}</option>}
          {tableOptions.map(o => <option key={o.key} value={o.key}>{o.label}</option>)}
        </select>
        <label style={{ fontWeight: 600, fontSize: '.85rem' }}>Tag by:</label>
        <select
          value={tag}
          onChange={e => setTag(e.target.value as 'rank' | 'category')}
          style={{ padding: '4px 8px', borderRadius: 4, border: '1px solid var(--color-border)' }}
        >
          <option value="category">Category</option>
          <option value="rank">Rank</option>
        </select>
        {byCategory ? (
          <select
            value={selectedCatSet}
            onChange={e => setSelectedCatSet(e.target.value)}
            style={{ padding: '4px 8px', borderRadius: 4, border: '1px solid var(--color-border)' }}
          >
            {categorySets.map(cs => (
              <option key={cs.name} value={cs.name}>{cs.label}</option>
            ))}
          </select>
        ) : (
          <select
            value={rank ?? ''}
            onChange={e => setRank(e.target.value)}
            style={{ padding: '4px 8px', borderRadius: 4, border: '1px solid var(--color-border)' }}
          >
            {ranks.map(r => <option key={r} value={r}>{r}</option>)}
          </select>
        )}
        {(subgroups ?? []).length >= 2 && (
          <>
            <label style={{ fontWeight: 600, fontSize: '.85rem' }}>Sub-group:</label>
            <select
              value={subgroup ?? ''}
              onChange={e => setSubgroup(e.target.value || null)}
              style={{ padding: '4px 8px', borderRadius: 4, border: '1px solid var(--color-border)' }}
            >
              <option value="">All</option>
              {(subgroups ?? []).map(sg => (
                <option key={sg} value={sg}>{sg}</option>
              ))}
            </select>
          </>
        )}
      </div>

      {/* Category summary table */}
      {summary && (
        <div style={{ marginBottom: 12, overflowX: 'auto' }}>
          <table style={{ borderCollapse: 'collapse', fontSize: '.82rem', width: '100%', maxWidth: 600 }}>
            <thead>
              <tr style={{ borderBottom: '2px solid var(--color-border)' }}>
                {byCategory && <th style={{ textAlign: 'center', padding: '4px 10px' }}>Colour</th>}
                <th style={{ textAlign: 'left', padding: '4px 10px' }}>{byCategory ? 'Category' : rank}</th>
                <th style={{ textAlign: 'right', padding: '4px 10px' }}>Rows</th>
                <th style={{ textAlign: 'right', padding: '4px 10px' }}>Reads</th>
                <th style={{ textAlign: 'right', padding: '4px 10px' }}>%</th>
              </tr>
            </thead>
            <tbody>
              {Object.entries(summary.categories).map(([name, stats]) => (
                <tr key={name} style={{ borderBottom: '1px solid var(--color-border)' }}>
                  {byCategory && <td style={{ padding: '3px 10px', textAlign: 'center' }}>
                    {savableCategories.has(name) ? (
                      <input
                        type="color"
                        value={effectiveColourMap[name] ?? UNASSIGNED_COLOUR}
                        onChange={e => setCategoryColour(name, e.target.value)}
                        title={`Colour for ${name}`}
                        disabled={catSetBusy}
                        style={{ width: 26, height: 20, padding: 0, border: '1px solid var(--color-border)',
                                 borderRadius: 4, background: 'none', cursor: 'pointer', verticalAlign: 'middle' }}
                      />
                    ) : (
                      // Catch-all bucket with a fixed swatch.
                      <span
                        title={`${name} (fixed colour)`}
                        style={{ display: 'inline-block', width: 26, height: 20,
                                 border: '1px solid var(--color-border)', borderRadius: 4,
                                 background: effectiveColourMap[name] ?? UNASSIGNED_COLOUR, verticalAlign: 'middle' }}
                      />
                    )}
                  </td>}
                  <td style={{ padding: '3px 10px', fontWeight: 600, color: 'var(--color-fg)' }}>
                    {name}
                  </td>
                  <td style={{ padding: '3px 10px', textAlign: 'right' }}>{stats.rows.toLocaleString()}</td>
                  <td style={{ padding: '3px 10px', textAlign: 'right' }}>{stats.reads.toLocaleString()}</td>
                  <td style={{ padding: '3px 10px', textAlign: 'right' }}>{stats.reads_percent.toFixed(2)}%</td>
                </tr>
              ))}
              <tr style={{ borderTop: '2px solid var(--color-border)', fontWeight: 600 }}>
                {byCategory && <td style={{ padding: '3px 10px' }}></td>}
                <td style={{ padding: '3px 10px' }}>Total</td>
                <td style={{ padding: '3px 10px', textAlign: 'right' }}>{summary.total_rows.toLocaleString()}</td>
                <td style={{ padding: '3px 10px', textAlign: 'right' }}>{summary.total_reads.toLocaleString()}</td>
                <td style={{ padding: '3px 10px', textAlign: 'right' }}>100%</td>
              </tr>
            </tbody>
          </table>

          {/* Colour action row: persist the pending overrides into the set */}
          {byCategory && <div style={{ display: 'flex', gap: 8, alignItems: 'center', marginTop: 10, flexWrap: 'wrap' }}>
            <button
              className="btn btn-primary"
              onClick={handleSaveCatSet}
              disabled={catSetBusy || !hasPendingColours}
              title="Overwrite the selected category set with the edited colours"
            >
              {catSetBusy ? 'Saving…' : 'Save'}
            </button>
            <button
              className="btn"
              onClick={() => setShowSaveAs(true)}
              disabled={catSetBusy}
              title="Save the edited colours as a new category set"
            >
              Save as…
            </button>
            {selectedCatSet !== 'default' && (
              <button
                className="btn btn-danger"
                onClick={handleDeleteCatSet}
                disabled={catSetBusy}
                title="Delete the selected category set"
              >
                Delete
              </button>
            )}
            {hasPendingColours && (
              <span style={{ fontSize: '.78rem', color: 'var(--color-muted-fg)' }}>Unsaved colour changes</span>
            )}
          </div>}
        </div>
      )}

      {/* Data table */}
      <DataTable
        key={`comp:${study}/${run}/${group ?? ''}/${effectiveSource}/${table}/${selectedCatSet}`}
        storageKey={`comp:${study}/${run}/${group ?? ''}/${effectiveSource}/${table}/${selectedCatSet}`}
        fetcher={fetcher}
        distinctFetcher={distinctFetcher}
      />

      {showSaveAs && (
        <NameDialog
          title="Save category set as"
          placeholder="New category-set name"
          onConfirm={handleSaveAsCatSet}
          onClose={() => setShowSaveAs(false)}
        />
      )}
    </div>
  )
}
