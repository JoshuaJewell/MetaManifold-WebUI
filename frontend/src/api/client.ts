import type {
  Study, StudySummary, Run,
  Job, JobStatus,
  TableMeta, TablePage, TableQuery, ColFilter, DistinctInfo,
  FilterPreset, ApplyPresetResult,
  ConfigMap,
  DatabaseEntry,
  ApiError,
  AnnotationSource,
  AnalysisRequest, ComparisonRequest, PermanovaResult,
  PublicationTable, PublicationTableRequest,
  TableDisplay,
  ChartRequest, CrossRunChartRequest, FacetChartRequest,
  CategorySet, CategorySetSaveRequest, CompositionBuildResult,
  CompositionSummaryRequest,
  VennRequest, VennResult,
  CompositionFilter, CompositionSet, CompositionLibraryDoc,
  PrimerDocument, PrimerSaveResult,
  DatabaseDocument, DatabaseSaveResult,
  TreeFileInfo, TreeFile, RunRef, ReadFunnelData, ReportItem, ReportKind,
  PlacementDoc, PlacementSummary, PlacementDetail, PlacementQueries, PhyloStep,
  ReferenceTreeDoc, ReferenceTreeSummary, ReferenceTreeDetail, AlignmentQC, TrimSettings,
} from './types'

import type { FigureDoc } from '../figure/types'
import { saveBlob } from '../utils/download'

// Base URL for the backend API. Empty string means same-origin.
// Set "apiBase" in config.json (e.g. "https://bioserver:8080") for split deployments.
let _apiBase = ''

/** Called once at startup from main.tsx to load runtime config. */
export async function loadConfig(): Promise<void> {
  try {
    const res = await fetch('/config.json')
    if (res.ok) {
      const cfg = await res.json()
      _apiBase = (cfg.apiBase as string ?? '').replace(/\/+$/, '')
    }
  } catch {
    // Missing or malformed config.json - default to same-origin
  }
}

/** Prepend the API base to a path that starts with "/" */
export function apiUrl(path: string): string {
  return _apiBase ? `${_apiBase}${path}` : path
}

async function rawPost(path: string, body: unknown): Promise<Response> {
  const res = await fetch(apiUrl(path), {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  })
  if (!res.ok) {
    const err: ApiError = await res.json().catch(() => ({
      error: 'network_error', message: res.statusText,
    }))
    throw Object.assign(new Error(err.message), { apiError: err, status: res.status })
  }
  return res
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const res = await fetch(apiUrl(path), {
    cache: 'no-store',
    headers: { 'Content-Type': 'application/json', ...init?.headers },
    ...init,
  })
  if (!res.ok) {
    const err: ApiError = await res.json().catch(() => ({
      error: 'network_error', message: res.statusText,
    }))
    throw Object.assign(new Error(err.message), { apiError: err, status: res.status })
  }
  if (res.status === 204) return undefined as T
  return res.json() as Promise<T>
}

const get  = <T>(path: string)                       => request<T>(path)
const post = <T>(path: string, body?: unknown)       => request<T>(path, { method: 'POST',  body: body !== undefined ? JSON.stringify(body) : undefined })
const patch = <T>(path: string, body?: unknown)      => request<T>(path, { method: 'PATCH', body: body !== undefined ? JSON.stringify(body) : undefined })
const put  = <T>(path: string, body?: unknown)       => request<T>(path, { method: 'PUT',  body: body !== undefined ? JSON.stringify(body) : undefined })
const del  = <T>(path: string)                       => request<T>(path, { method: 'DELETE' })

/** PUT a text body; the response is JSON. */
async function putText<T>(path: string, text: string): Promise<T> {
  const res = await fetch(apiUrl(path), { method: 'PUT', headers: { 'Content-Type': 'text/plain' }, body: text })
  if (!res.ok) {
    const err: ApiError = await res.json().catch(() => ({ error: 'network_error', message: res.statusText }))
    throw Object.assign(new Error(err.message), { apiError: err, status: res.status })
  }
  return res.json() as Promise<T>
}

/** GET a text body, or null when there is none yet. */
async function getText(path: string): Promise<string | null> {
  const res = await fetch(apiUrl(path), { cache: 'no-store' })
  return res.ok ? res.text() : null
}

/** Append ?group=X query parameter when group is provided. */
const gq = (group?: string | null) =>
  group ? `?group=${encodeURIComponent(group)}` : ''

export const api = {
  studies: {
    list:   ()               => get<StudySummary[]>('/api/v1/studies'),
    get:    (study: string)  => get<Study>(`/api/v1/studies/${study}`),
    create: (name: string)   => post<Study>('/api/v1/studies', { name }),
    rename: (study: string, name: string) => post<Study>(`/api/v1/studies/${study}/rename`, { name }),
    delete: (study: string)  => del<{ deleted: string }>(`/api/v1/studies/${study}`),
  },

  groups: {
    create: (study: string, name: string)               => post<{ study: string; name: string }>(`/api/v1/studies/${study}/groups`, { name }),
    rename: (study: string, group: string, name: string) => post<{ study: string; name: string }>(`/api/v1/studies/${study}/groups/${group}/rename`, { name }),
    delete: (study: string, group: string)              => del<{ deleted: string }>(`/api/v1/studies/${study}/groups/${group}`),
  },

  runs: {
    list:      (study: string)                          => get<Run[]>(`/api/v1/studies/${study}/runs`),
    get:       (study: string, run: string, group?: string | null) => get<Run>(`/api/v1/studies/${study}/runs/${run}${gq(group)}`),
    listGroup: (study: string, group: string)            => get<Run[]>(`/api/v1/studies/${study}/groups/${group}/runs`),
    create:    (study: string, name: string, group?: string | null) =>
      post<RunRef>(`/api/v1/studies/${study}/runs`, { name, ...(group ? { group } : {}) }),
    rename:    (study: string, run: string, name: string, group?: string | null) => post<RunRef>(`/api/v1/studies/${study}/runs/${run}/rename${gq(group)}`, { name }),
    delete:    (study: string, run: string, group?: string | null) => del<{ deleted: string }>(`/api/v1/studies/${study}/runs/${run}${gq(group)}`),
  },

  pipeline: {
    runStudy: (study: string)                       => post<Job>(`/api/v1/studies/${study}/pipeline`),
    runRun:   (study: string, run: string, group?: string | null) => post<Job>(`/api/v1/studies/${study}/runs/${run}/pipeline${gq(group)}`),
    runStage: (study: string, run: string, stage: string, group?: string | null) =>
                                                       post<Job>(`/api/v1/studies/${study}/runs/${run}/stages/${stage}${gq(group)}`),
  },

  jobs: {
    list:   (opts?: { study?: string; status?: JobStatus }) => {
      const p = new URLSearchParams()
      if (opts?.study)  p.set('study',  opts.study)
      if (opts?.status) p.set('status', opts.status)
      const qs = p.size ? `?${p}` : ''
      return get<Job[]>(`/api/v1/jobs${qs}`)
    },
    cancel: (id: string) => del<void>(`/api/v1/jobs/${id}`),
  },

  results: {
    runTables:    (study: string, run: string, group?: string | null) => get<TableMeta[]>(`/api/v1/studies/${study}/runs/${run}/results/tables${gq(group)}`),
    runTable:     (study: string, run: string, id: string, q: TableQuery, group?: string | null) =>
      post<TablePage>(`/api/v1/studies/${study}/runs/${run}/results/tables/${id}/query${gq(group)}`, q),
    distinctValues: (study: string, run: string, id: string, column: string, activeFilters?: Record<string, ColFilter>, group?: string | null, keywordFilter?: string) =>
      post<DistinctInfo>(`/api/v1/studies/${study}/runs/${run}/results/tables/${id}/distinct/${column}${gq(group)}`,
        {
          ...(activeFilters ? { colFilters: activeFilters } : {}),
          ...(keywordFilter ? { filter: keywordFilter } : {}),
        }),
    saveTable: (study: string, run: string, id: string, name: string,
                colFilters?: Record<string, ColFilter>, sortBy?: string, sortDir?: string, group?: string | null) =>
      post<{ name: string; path: string; rows: number }>(
        `/api/v1/studies/${study}/runs/${run}/results/tables/${id}/save${gq(group)}`,
        { name, colFilters, sortBy, sortDir }),
    deleteTable: (study: string, run: string, id: string, group?: string | null) =>
      del<{ deleted: string }>(`/api/v1/studies/${study}/runs/${run}/results/tables/${id}${gq(group)}`),
    otuMembers: (study: string, run: string, otu: string, group?: string | null) =>
      get<{ otu: string; columns: string[]; rows: Record<string, unknown>[] }>(
        `/api/v1/studies/${study}/runs/${run}/results/otu-members/${otu}${gq(group)}`),
    otuCounts: (study: string, run: string, group?: string | null) =>
      get<{ counts: Record<string, number> }>(
        `/api/v1/studies/${study}/runs/${run}/results/otu-counts${gq(group)}`),
    qcOutputs: (study: string, run: string, group?: string | null) =>
      get<{ has_report: boolean; report_url: string | null }>(
        `/api/v1/studies/${study}/runs/${run}/results/qc${gq(group)}`),
    dada2Outputs: (study: string, run: string, group?: string | null) =>
      get<{ figures: { name: string; label: string; url: string }[]; has_stats: boolean; logs: { name: string; url: string }[]; config: Record<string, unknown> }>(
        `/api/v1/studies/${study}/runs/${run}/results/dada2${gq(group)}`),
    dada2Stats: (study: string, run: string, group?: string | null) =>
      get<{ columns: string[]; rows: Record<string, unknown>[] }>(
        `/api/v1/studies/${study}/runs/${run}/results/dada2/stats${gq(group)}`),
    exportTable: async (study: string, run: string, id: string,
                        colFilters?: Record<string, ColFilter>, sortBy?: string, sortDir?: string, group?: string | null,
                        display?: TableDisplay) => {
      const res = await rawPost(
        `/api/v1/studies/${study}/runs/${run}/results/tables/${id}/export${gq(group)}`,
        { colFilters, sortBy, sortDir, filename: `${study}_${run}_filtered.xlsx`, ...display })
      saveBlob(await res.blob(), `${study}_${run}_filtered.xlsx`)
    },
  },

  presets: {
    list: () => get<FilterPreset[]>('/api/v1/filter-presets'),
    apply: (study: string, run: string, table: string, preset: string, group?: string | null) =>
      post<ApplyPresetResult>(`/api/v1/studies/${study}/runs/${run}/results/tables/${table}/apply-preset${gq(group)}`, { preset }),
    save: (name: string, filters: Record<string, ColFilter>, description?: string) =>
      post<FilterPreset>(`/api/v1/filter-presets/${name}`, { filters, description }),
    delete: (name: string) => del<{ deleted: string }>(`/api/v1/filter-presets/${name}`),
  },

  primers: {
    list:     () => get<string[]>('/api/v1/primers'),
    document: () => get<PrimerDocument>('/api/v1/primers/document'),
    save:     (doc: PrimerDocument) => put<PrimerSaveResult>('/api/v1/primers', doc),
  },

  config: {
    getDefault:    ()                                => get<ConfigMap>('/api/v1/config'),
    patchDefault:  (body: Record<string, unknown>)   => patch<ConfigMap>('/api/v1/config', body),
    deleteDefault: (key: string)                     => del<ConfigMap>(`/api/v1/config/${encodeURIComponent(key)}`),
    getStudy:    (study: string)                     => get<ConfigMap>(`/api/v1/studies/${study}/config`),
    patchStudy:  (study: string, body: Record<string, unknown>) => patch<ConfigMap>(`/api/v1/studies/${study}/config`, body),
    deleteStudy: (study: string, key: string)        => del<ConfigMap>(`/api/v1/studies/${study}/config/${key}`),
    getRun:      (study: string, run: string, group?: string | null) => get<ConfigMap>(`/api/v1/studies/${study}/runs/${run}/config${gq(group)}`),
    patchRun:    (study: string, run: string, body: Record<string, unknown>, group?: string | null) =>
                                                        patch<ConfigMap>(`/api/v1/studies/${study}/runs/${run}/config${gq(group)}`, body),
    deleteRun:   (study: string, run: string, key: string, group?: string | null) =>
                                                        del<ConfigMap>(`/api/v1/studies/${study}/runs/${run}/config/${key}${gq(group)}`),
    studyOverrides: (study: string) => get<Record<string, string[]>>(`/api/v1/studies/${study}/config/overrides`),
    groupOverrides: (study: string, group: string) => get<Record<string, string[]>>(`/api/v1/studies/${study}/groups/${group}/config/overrides`),
    getGroup:    (study: string, group: string)          => get<ConfigMap>(`/api/v1/studies/${study}/groups/${group}/config`),
    patchGroup:  (study: string, group: string, body: Record<string, unknown>) =>
                                                        patch<ConfigMap>(`/api/v1/studies/${study}/groups/${group}/config`, body),
    deleteGroup: (study: string, group: string, key: string) =>
                                                        del<ConfigMap>(`/api/v1/studies/${study}/groups/${group}/config/${key}`),
  },

  analysis: {
    alpha:         (study: string, run: string, body: AnalysisRequest, group?: string | null) =>
                     post<unknown>(`/api/v1/studies/${study}/runs/${run}/analysis/alpha${gq(group)}`, body),
    chart:         (study: string, run: string, body: ChartRequest, group?: string | null) =>
                     post<unknown>(`/api/v1/studies/${study}/runs/${run}/analysis/chart${gq(group)}`, body),
    readFunnel:    (study: string, run: string, group?: string | null) =>
                     get<ReadFunnelData>(`/api/v1/studies/${study}/runs/${run}/analysis/read-funnel${gq(group)}`),
    ranks:         (study: string, run: string, opts?: { table?: string; group?: string | null; source?: AnnotationSource }) => {
                     const query = new URLSearchParams()
                     if (opts?.group) query.set('group', opts.group)
                     if (opts?.table) query.set('table', opts.table)
                     if (opts?.source) query.set('source', opts.source)
                     const qs = query.size ? `?${query.toString()}` : ''
                     return get<string[]>(`/api/v1/studies/${study}/runs/${run}/analysis/ranks${qs}`)
                   },
    compareAlpha:  (study: string, body: ComparisonRequest) =>
                     post<unknown>(`/api/v1/studies/${study}/analysis/alpha`, body),
    chartCompare:  (study: string, body: CrossRunChartRequest) =>
                     post<unknown>(`/api/v1/studies/${study}/analysis/chart`, body),
    chartFacet:    (study: string, body: FacetChartRequest) =>
                     post<unknown>(`/api/v1/studies/${study}/analysis/chart-facet`, body),
    nmds:          (study: string, body: ComparisonRequest) =>
                     post<unknown>(`/api/v1/studies/${study}/analysis/nmds`, body),
    permanova:     (study: string, body: ComparisonRequest) =>
                     post<PermanovaResult>(`/api/v1/studies/${study}/analysis/permanova`, body),
    capabilities:  () => get<{ r_available: boolean }>('/api/v1/capabilities'),
    venn:          (study: string, body: VennRequest) =>
                     post<VennResult>(`/api/v1/studies/${study}/analysis/venn`, body),
    publicationTable: (study: string, body: PublicationTableRequest) =>
                     post<{ table: PublicationTable }>(`/api/v1/studies/${study}/analysis/publication-tables`, body),
    /** The formatted .xlsx of the table, as a blob. */
    publicationTableBlob: async (study: string, body: PublicationTableRequest) =>
      (await rawPost(`/api/v1/studies/${study}/analysis/publication-tables`, { ...body, format: 'xlsx' })).blob(),
    /** Download the table as a formatted .xlsx or as CSV. */
    downloadPublicationTable: async (study: string, body: PublicationTableRequest, format: 'xlsx' | 'csv') => {
      const res = await rawPost(`/api/v1/studies/${study}/analysis/publication-tables`, { ...body, format })
      const name = /filename="([^"]+)"/.exec(res.headers.get('Content-Disposition') ?? '')?.[1]
      saveBlob(await res.blob(), name ?? `${study}_table.${format}`)
    },
  },

  composition: {
    categorySets: () => get<CategorySet[]>('/api/v1/category-sets'),
    saveCategorySet: (name: string, body: CategorySetSaveRequest) =>
      post<CategorySet>(`/api/v1/category-sets/${encodeURIComponent(name)}`, body),
    deleteCategorySet: (name: string) =>
      del<{ deleted: string }>(`/api/v1/category-sets/${encodeURIComponent(name)}`),
    summary: (study: string, run: string, body: CompositionSummaryRequest, group?: string | null) =>
      post<CompositionBuildResult>(
        `/api/v1/studies/${study}/runs/${run}/composition/summary${gq(group)}`, body),
    query: (study: string, run: string, source: AnnotationSource, q: TableQuery,
            group?: string | null, categorySet?: string, table?: string) =>
      post<TablePage>(
        `/api/v1/studies/${study}/runs/${run}/composition/${source}/query${gq(group)}`,
        { ...q, ...(categorySet ? { category_set: categorySet } : {}), ...(table ? { table } : {}) }),
    distinct: (study: string, run: string, source: AnnotationSource, column: string,
               activeFilters?: Record<string, ColFilter>, group?: string | null,
               categorySet?: string, keywordFilter?: string, table?: string) =>
      post<DistinctInfo>(
        `/api/v1/studies/${study}/runs/${run}/composition/${source}/distinct/${column}${gq(group)}`,
        {
          ...(activeFilters ? { colFilters: activeFilters } : {}),
          ...(categorySet ? { category_set: categorySet } : {}),
          ...(keywordFilter ? { filter: keywordFilter } : {}),
          ...(table ? { table } : {}),
        }),
    library:      () => get<CompositionLibraryDoc>('/api/v1/composition'),
    saveFilter:   (name: string, body: CompositionFilter) =>
      post<CompositionFilter>(`/api/v1/composition/filters/${encodeURIComponent(name)}`, body),
    deleteFilter: (name: string) =>
      del<{ deleted: string }>(`/api/v1/composition/filters/${encodeURIComponent(name)}`),
    saveSet:      (name: string, body: CompositionSet) =>
      post<CompositionSet>(`/api/v1/composition/sets/${encodeURIComponent(name)}`, body),
    deleteSet:    (name: string) =>
      del<{ deleted: string }>(`/api/v1/composition/sets/${encodeURIComponent(name)}`),
  },


  databases: {
    list:     ()             => get<DatabaseEntry[]>('/api/v1/databases'),
    download: (key: string)  => post<Job>(`/api/v1/databases/${key}/download`),
    document: () => get<DatabaseDocument>('/api/v1/databases/document'),
    save:     (doc: DatabaseDocument) => put<DatabaseSaveResult>('/api/v1/databases', doc),
  },

  report: {
    list:    (study: string) => get<ReportItem[]>(`/api/v1/studies/${study}/report`),
    /** Upload a file as a new item; the body is sent as raw bytes. */
    add:     async (study: string, kind: ReportKind, title: string, ext: string, blob: Blob) => {
      const q = new URLSearchParams({ kind, title, ext })
      const res = await fetch(apiUrl(`/api/v1/studies/${study}/report?${q}`), {
        method: 'POST', headers: { 'Content-Type': 'application/octet-stream' }, body: blob,
      })
      if (!res.ok) {
        const err: ApiError = await res.json().catch(() => ({ error: 'network_error', message: res.statusText }))
        throw Object.assign(new Error(err.message), { apiError: err, status: res.status })
      }
      return res.json() as Promise<ReportItem>
    },
    rename:  (study: string, id: string, title: string) =>
      patch<ReportItem>(`/api/v1/studies/${study}/report/${id}`, { title }),
    reorder: (study: string, ids: string[]) =>
      put<ReportItem[]>(`/api/v1/studies/${study}/report/order`, { ids }),
    remove:  (study: string, id: string) => del<{ deleted: string }>(`/api/v1/studies/${study}/report/${id}`),
    fileUrl: (study: string, id: string) => apiUrl(`/api/v1/studies/${study}/report/${id}/file`),
    exportUrl: (study: string) => apiUrl(`/api/v1/studies/${study}/report/export`),
  },

  figures: {
    list:   (study: string) => get<{ id: string; title: string; modified: string }[]>(`/api/v1/studies/${study}/figures`),
    get:    (study: string, id: string) => get<FigureDoc>(`/api/v1/studies/${study}/figures/${id}`),
    create: (study: string, doc: Omit<FigureDoc, 'id'>) => post<FigureDoc>(`/api/v1/studies/${study}/figures`, doc),
    save:   (study: string, doc: FigureDoc) => put<FigureDoc>(`/api/v1/studies/${study}/figures/${doc.id}`, doc),
    remove: (study: string, id: string) => del<{ deleted: string }>(`/api/v1/studies/${study}/figures/${id}`),
    /** Save the rendered page beside the layout; the body is the PNG. */
    saveRaster: async (study: string, id: string, png: Blob) => {
      const res = await fetch(apiUrl(`/api/v1/studies/${study}/figures/${id}/raster`), {
        method: 'PUT', headers: { 'Content-Type': 'image/png' }, body: png,
      })
      if (!res.ok) {
        const err: ApiError = await res.json().catch(() => ({ error: 'network_error', message: res.statusText }))
        throw Object.assign(new Error(err.message), { apiError: err, status: res.status })
      }
      return res.json() as Promise<{ file: string }>
    },
  },

  referenceTrees: {
    list:    () => get<ReferenceTreeSummary[]>('/api/v1/reference-trees'),
    get:     (id: string) => get<ReferenceTreeDetail>(`/api/v1/reference-trees/${id}`),
    create:  (name: string) => post<ReferenceTreeDoc>('/api/v1/reference-trees', { name }),
    save:    (id: string, patch: Partial<Pick<ReferenceTreeDoc, 'name' | 'description' | 'overrides'>>) =>
      put<ReferenceTreeDoc>(`/api/v1/reference-trees/${id}`, patch),
    remove:  (id: string) => del<{ deleted: string }>(`/api/v1/reference-trees/${id}`),
    uploadFasta: (id: string, text: string) => putText<{ count: number; renamed: number }>(`/api/v1/reference-trees/${id}/fasta`, text),
    fastaUrl: (id: string) => apiUrl(`/api/v1/reference-trees/${id}/fasta`),
    treefileUrl: (id: string) => apiUrl(`/api/v1/reference-trees/${id}/treefile`),
    tree:    (id: string) => get<TreeFile>(`/api/v1/reference-trees/${id}/tree`),
    saveView: (id: string, view: unknown) => put<{ saved: boolean }>(`/api/v1/reference-trees/${id}/tree/view`, view),
    run:     (id: string) => post<Job>(`/api/v1/reference-trees/${id}/run`),
    log:     (id: string, step: PhyloStep) => getText(`/api/v1/reference-trees/${id}/log/${step}`),
    qc:      (id: string, step: PhyloStep) => get<AlignmentQC>(`/api/v1/reference-trees/${id}/qc/${step}`),
    alignment: (id: string, which: 'raw' | 'trimmed') => getText(`/api/v1/reference-trees/${id}/alignment/${which}`),
    trimPreview: (id: string, trim: TrimSettings) => post<AlignmentQC>(`/api/v1/reference-trees/${id}/trim-preview`, trim),
  },

  placements: {
    list:    (study: string) => get<PlacementSummary[]>(`/api/v1/studies/${study}/placements`),
    get:     (study: string, id: string) => get<PlacementDetail>(`/api/v1/studies/${study}/placements/${id}`),
    create:  (study: string, name: string, reference: string | null) =>
      post<PlacementDoc>(`/api/v1/studies/${study}/placements`, { name, reference }),
    save:    (study: string, id: string, patch: Partial<Pick<PlacementDoc, 'name' | 'reference' | 'queries' | 'overrides'>>) =>
      put<PlacementDoc>(`/api/v1/studies/${study}/placements/${id}`, patch),
    remove:  (study: string, id: string) => del<{ deleted: string }>(`/api/v1/studies/${study}/placements/${id}`),
    uploadQueries: (study: string, id: string, text: string) =>
      putText<{ count: number; renamed: number }>(`/api/v1/studies/${study}/placements/${id}/fasta/queries`, text),
    queriesUrl: (study: string, id: string) => apiUrl(`/api/v1/studies/${study}/placements/${id}/fasta/queries`),
    preview: (study: string, queries: PlacementQueries) =>
      post<{ count: number; per_run: { run: string; group: string | null; subgroups: string[]; count: number }[] }>(
        `/api/v1/studies/${study}/placements/preview`, queries),
    run:     (study: string, id: string) => post<Job>(`/api/v1/studies/${study}/placements/${id}/run`),
    log:     (study: string, id: string, step: PhyloStep) => getText(`/api/v1/studies/${study}/placements/${id}/log/${step}`),
    qc:      (study: string, id: string, step: PhyloStep) => get<AlignmentQC>(`/api/v1/studies/${study}/placements/${id}/qc/${step}`),
    alignment: (study: string, id: string, which: 'raw' | 'trimmed') =>
      getText(`/api/v1/studies/${study}/placements/${id}/alignment/${which}`),
    trimPreview: (study: string, id: string, trim: TrimSettings) =>
      post<AlignmentQC>(`/api/v1/studies/${study}/placements/${id}/trim-preview`, trim),
  },

  trees: {
    list:     (study: string) => get<TreeFileInfo[]>(`/api/v1/studies/${study}/trees`),
    get:      (study: string, file: string) =>
      get<TreeFile>(`/api/v1/studies/${study}/trees/${encodeURIComponent(file)}`),
    upload:   (study: string, file: string, content: string, overwrite = false) =>
      post<{ file: string }>(`/api/v1/studies/${study}/trees`, { file, content, overwrite }),
    saveView: (study: string, file: string, view: unknown) =>
      put<{ saved: boolean }>(`/api/v1/studies/${study}/trees/${encodeURIComponent(file)}/view`, view),
    delete:   (study: string, file: string) =>
      del<{ deleted: string }>(`/api/v1/studies/${study}/trees/${encodeURIComponent(file)}`),
  },
}
