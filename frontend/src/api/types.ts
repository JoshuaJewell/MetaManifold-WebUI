export interface StudySummary {
  name:             string
  run_count:        number
  group_count:      number
  active_job_count: number
}

export interface Study extends StudySummary {
  runs:   string[]
  groups: string[]
}

export type StageStatus = 'not_started' | 'running' | 'complete' | 'stale' | 'disabled'

export interface StageInfo {
  status:      StageStatus
  last_run:    string | null
  stale_keys?: string[]
}

export interface RunStages {
  fastqc:         StageInfo
  cutadapt:       StageInfo
  dada2_denoise:  StageInfo
  dada2_classify: StageInfo
  cdhit:          StageInfo
  swarm:          StageInfo
  vsearch:        StageInfo
  merge_taxa:     StageInfo
}

export interface Run {
  name:         string
  study:        string
  sample_count: number
  samples:      string[]
  stages:       RunStages
  pooled:       boolean
  subgroups:    string[]
  /** Classifiers whose taxonomy the run's merged table carries; null when the table could not be read. */
  taxonomy_sources: AnnotationSource[] | null
}

/** What the run create and rename routes return. */
export interface RunRef {
  study: string
  name:  string
  group?: string | null
}

export type JobStatus = 'queued' | 'running' | 'complete' | 'failed' | 'cancelled'
export type JobType   = 'pipeline' | 'stage' | 'db_download'

export interface Job {
  id:          string
  type:        JobType
  study:       string | null
  run:         string | null
  group?:      string | null
  stage:       string | null
  status:      JobStatus
  created_at:  string
  finished_at: string | null
  message:     string | null
}

export interface TableMeta {
  id:           string
  label:        string
  rows:         number
  total_reads:  number
  n_samples:    number
}

export interface ColFilter {
  text?:     string
  include?:  string[]
  exclude?:  string[]
  min?:      number
  max?:      number
  /** Only on the SAMPLE_READS_FILTER_KEY entry: which rows a sample's total is measured over. */
  basis?:    SampleReadsBasis
}

/**
 * Reserved colFilters key carrying the sample total-read bounds. A sample is a
 * column, so this entry drops whole sample columns (from the table, its read
 * totals, charts, save and export) rather than rows. Mirrors
 * SAMPLE_READS_FILTER_KEY in src/server/routes/duckdb_helpers.jl.
 */
export const SAMPLE_READS_FILTER_KEY = '__sample_reads__'

/** 'filtered': reads left after the row filters (default). 'raw': the whole library. */
export type SampleReadsBasis = 'filtered' | 'raw'

export interface TableQuery {
  page:        number
  perPage:     number
  filter?:     string
  sortBy?:     string
  sortDir?:    'asc' | 'desc'
  colFilters?: Record<string, ColFilter>
}

export interface DistinctText {
  column: string
  type:   'text'
  values: string[]
  count:  number
}

export interface DistinctNumeric {
  column: string
  type:   'numeric'
  min:    number
  max:    number
  count:  number
  sum:    number
  mean:   number
  median: number
  q1:     number
  q3:     number
}

export type DistinctInfo = DistinctText | DistinctNumeric

export interface TablePage {
  total:                  number
  total_unfiltered:       number
  total_reads:            number
  total_reads_unfiltered: number
  page:                   number
  per_page:               number
  columns:                string[]
  sample_count_columns:   string[]
  /** Sample columns dropped by the sample read-count filter. */
  excluded_samples?:      string[]
  /** Largest value of each count column over the filtered rows. */
  count_max?:             Record<string, number>
  rows:                   Record<string, unknown>[]
}

/** Where a setting's value comes from. placement and tree are a placement's or reference tree's own overrides. */
export type ConfigSource = 'default' | 'study' | 'group' | 'run' | 'placement' | 'tree'

export type ConfigMap = Record<string, { value: unknown; source: ConfigSource }>

export interface DatabaseEntry {
  key:              string
  label:            string
  dada2_available:  boolean
  vsearch_available: boolean
}

export interface FilterPreset {
  name:        string
  label:       string
  file:        string
  description: string
}

export interface ApplyPresetResult {
  preset:      string
  rows_before: number
  rows_after:  number
  filters:     Record<string, ColFilter>
}

export interface AnalysisRequest {
  table: string
  source?: AnnotationSource
  colFilters?: Record<string, ColFilter>
  prefix?: string | null
}

export interface ChartRequest {
  table:        string
  /** Classifier to read ranks and categories from; the configured one when absent. */
  source?:      AnnotationSource
  tag:          'rank' | 'category'
  value:        string
  relative?:    boolean
  mode?:        'stacked' | 'grouped'
  subgroup?:    string | null
  top_n?:       number
  /** Show samples with no reads as blank slots on the axis. */
  keep_empty?:  boolean
  colFilters?:  Record<string, ColFilter>
}

export interface CrossRunChartRequest extends ChartRequest {
  runs: ComparisonRunSpec[]
}

/** The dimensions a facet axis can be keyed on. */
export const FACET_DIMENSIONS = ['run', 'group', 'subgroup'] as const
export type FacetDimension = (typeof FACET_DIMENSIONS)[number]

/** A grid of composition charts: one panel per (rows value, cols value) pair. */
export interface FacetChartRequest extends Omit<ChartRequest, 'subgroup'> {
  runs: ComparisonRunSpec[]
  rows: FacetDimension
  cols: FacetDimension
}

export interface ComparisonRunSpec {
  run: string
  group?: string | null
  prefix?: string | null
  source?: AnnotationSource
}

/** Expand pooled runs into per-subgroup ComparisonRunSpecs. */
export function expandRunSpecs(
  items: Pick<Run, 'name' | 'pooled' | 'subgroups'>[],
  group?: string | null,
): ComparisonRunSpec[] {
  return items.flatMap(item =>
    item.pooled && item.subgroups.length > 0
      ? item.subgroups.map(prefix => ({ run: item.name, group: group ?? null, prefix }))
      : [{ run: item.name, group: group ?? null }]
  )
}

/** Whether a run's tables carry taxonomy from `source`. Unknown counts as yes. */
export function carriesSource(run: Pick<Run, 'taxonomy_sources'>, source: AnnotationSource): boolean {
  return run.taxonomy_sources == null || run.taxonomy_sources.includes(source)
}

/** One spec per distinct run and group. */
export function uniqueRuns(specs: ComparisonRunSpec[]): ComparisonRunSpec[] {
  const seen = new Set<string>()
  return specs.filter(s => {
    const key = `${s.group ?? ''}/${s.run}`
    if (seen.has(key)) return false
    seen.add(key)
    return true
  })
}

export interface ComparisonRequest extends AnalysisRequest {
  runs: ComparisonRunSpec[]
  aggregate?: boolean
}

export interface PermanovaResult {
  text: string
  r2: number | null
  f_statistic: number | null
  p_value: number | null
  /** True when permutations were restricted within individuals. */
  blocked?: boolean
  /** One row per model term, tested sequentially in formula order. */
  terms?: { term: string; r2: number; f_statistic: number | null; p_value: number | null }[]
  /** PERMDISP on the same distances: do the groups differ in spread? */
  dispersion?: PermdispResult | { error: string }
}

export interface PermdispResult {
  groups: { group: string; mean: number; median: number }[]
  f_statistic: number
  df: [number, number]
  p_value: number
  text: string
}

export type PublicationCellKind = 'label' | 'int' | 'pct' | 'pp'

export interface PublicationTable {
  title: string
  /** Spanning group headings above the column labels, outermost first. */
  header_rows: { label: string; span: number }[][]
  columns: { label: string; kind: PublicationCellKind }[]
  rows: (string | number | null)[][]
  /** Totals, set off from the body by a rule. */
  footer: (string | number | null)[][]
  notes: string[]
  /** Heatmap background per body cell (hex), when a heatmap was requested. */
  fills?: (string | null)[][]
}

/** `column` grades each column on its own; `table` shares one scale per kind. */
export type HeatmapMode = 'none' | 'column' | 'table'
export type PublicationHeatmap = HeatmapMode

/** Count-cell display options for a results table, also applied on export. */
export interface TableDisplay {
  heatmap:    HeatmapMode
  hide_zeros: boolean
}

export type PublicationValue = 'asvs' | 'reads' | 'pct'

export interface PublicationTableRequest {
  runs: ComparisonRunSpec[]
  /** Results table the counts come from. */
  table: string
  rows: 'rank' | 'category'
  rank: string
  /** With rows 'category': the set that sorts ASVs into categories. */
  category_set: string
  columns: 'run' | 'subgroup' | 'sample'
  values: PublicationValue[]
  /** A percentage-point difference for runs with exactly two sub-groups. */
  difference: boolean
  /** Empty uses a title describing the table. */
  title: string
  heatmap?: PublicationHeatmap
  /** Blank exact zeros in the table body (not totals). */
  hide_zeros?: boolean
}

export interface ApiError {
  error:   string
  message: string
  detail?: string
}

export type AnnotationSource = 'VSEARCH' | 'DADA2'

export interface CategoryInfo {
  name:    string
  colour?: string
}

export interface CategorySet {
  name:               string
  label:              string
  description:        string
  categories:         CategoryInfo[]
  unassigned_colour?: string
}

export interface CompositionCategoryStats {
  rows:          number
  reads:         number
  reads_percent: number
}

export interface CompositionBuildResult {
  table:        string
  source:       AnnotationSource
  category_set: string
  tag?:         'rank' | 'category'
  /** Category set or rank the labels come from. */
  value?:       string
  total_rows:   number
  total_reads:  number
  categories:   Record<string, CompositionCategoryStats>
}

export interface VennSet {
  name: string
  taxa: string[]
}

export interface VennResult {
  sets: VennSet[]
  rank: string
}

export interface VennRequest {
  runs: ComparisonRunSpec[]
  table: string
  rank: string
  source?: AnnotationSource
}

/** Exactly two conditions; the first is the reference. */
export interface DifferentialRequest {
  runs: ComparisonRunSpec[]
  table: string
  rank: string
  aggregate?: boolean
}

/** `ok` and `boundary` rows carry statistics; `failed` and `filtered` rows say why not. */
export type DifferentialStatus = 'ok' | 'boundary' | 'failed' | 'filtered'

export interface DifferentialRow {
  taxon: string
  status: DifferentialStatus
  note: string
  /**
   * nb_glm: natural-log fold change, contrast over reference. clr_welch: mean
   * CLR of contrast minus mean CLR of reference, which is not a fold change.
   */
  estimate: number | null
  /** nb_glm only; null for clr_welch. */
  log2_fold_change: number | null
  standard_error: number | null
  /** z (nb_glm) or Welch's t (clr_welch). */
  statistic: number | null
  pvalue: number | null
  /** Benjamini-Hochberg adjusted; null outside the tested family. */
  padj: number | null
  /** nb_glm only. */
  dispersion_theta: number | null
  prevalence: number
}

export type DifferentialMethod = 'nb_glm' | 'clr_welch'

/** Which row field holds the effect size, and what to call it. */
export interface DifferentialEffect {
  key: 'log2_fold_change' | 'estimate'
  label: string
}

export interface DifferentialResult {
  status: 'ok' | 'partial'
  method: string
  rank: string
  table: string
  groups: { reference: string; contrast: string }
  n_samples: Record<string, number>
  effect: DifferentialEffect
  /** Empty for clr_welch, which takes no library-size offset. */
  size_factors: { sample: string; group: string; factor: number }[]
  config: {
    method: DifferentialMethod
    min_prevalence: number
    /** nb_glm only. */
    offset?: 'tss' | 'rle'
    /** clr_welch only. */
    replacement_delta?: number
  }
  diagnostics: {
    n_taxa: number; n_tested: number; n_failed: number; n_boundary: number; n_filtered: number
    /** clr_welch only. */
    zero_replacement?: {
      method: string; delta: number; zeros_replaced: number; n_taxa_in_composition: number
      n_taxa_unobserved: number; median_imputed_fraction: number; max_imputed_fraction: number
      warnings: string[]
    }
  }
  rows: DifferentialRow[]
  figure: unknown
}

export interface CategorySetSaveRequest {
  base:         string
  colours:      Record<string, string>
  label?:       string
  description?: string
}

export interface CompositionSummaryRequest {
  category_set: string
  source?:      AnnotationSource
  subgroup?:    string | null
  /** Results table to summarise; the backend defaults to merged. */
  table?:       string
  /** Break down by category set (default) or rank, as in ChartRequest. */
  tag?:         'rank' | 'category'
  value?:       string
}

export interface CompositionFilterRule {
  column?:  string
  type?:    'include' | 'min' | 'max'
  values?:  string[]
  value?:   number
  // Pattern rules are preserved on round-trip but not editable in the first-pass UI.
  pattern?: string
  action?:  string
  regex?:   boolean
}

export interface CompositionFilter {
  databases?:    string[]
  remove_empty?: string[]
  filters?:      CompositionFilterRule[]
}

export interface CompositionCategory extends CategoryInfo {
  // An absent filter marks the catch-all category.
  filter?: string
}

export interface CompositionSet {
  label?:             string
  description?:       string
  unassigned_colour?: string
  categories:         CompositionCategory[]
}

export interface CompositionLibraryDoc {
  filters: Record<string, CompositionFilter>
  sets:    Record<string, CompositionSet>
}

export interface PrimerPair {
  name:    string
  forward: string
  reverse: string
}

export interface PrimerDocument {
  Forward: Record<string, string>
  Reverse: Record<string, string>
  Pairs:   PrimerPair[]
}

export interface PrimerWarning {
  pair:          string
  referenced_by: string[]
}

export interface PrimerSaveResult {
  document: PrimerDocument
  warnings: PrimerWarning[]
}

export interface DatabaseFormat {
  uri:          string
  local:        string | null
  // dada2 only: a pre-existing path on the remote taxonomy host.
  remote_path?: string | null
}

export interface DatabaseCorrectionValue {
  from: string
  to:   string
}

export interface DatabaseCorrection {
  source: string
  target: string
  // Ordered rows: a map keyed by edited text would collapse two rows as soon
  // as one key is typed into the other.
  values: DatabaseCorrectionValue[]
}

export interface DatabaseDocumentEntry {
  key:            string
  label:          string
  dada2:          DatabaseFormat
  vsearch:        DatabaseFormat
  levels:         string[]
  vsearch_format: string
  corrections:    DatabaseCorrection[]
}

export interface DatabaseDocument {
  dir:       string
  databases: DatabaseDocumentEntry[]
}

export interface DatabaseWarning {
  kind:             'database_removed' | 'levels_changed' | 'release_mismatch'
  database:         string
  used_by?:         string[]
  dada2_version?:   string
  vsearch_version?: string
}

export interface DatabaseSaveResult {
  document: DatabaseDocument
  warnings: DatabaseWarning[]
}

export interface TreeFileInfo {
  file:     string
  format:   'jplace' | 'newick'
  size:     number
  modified: string
  has_view: boolean
}

export interface TreeFile {
  file:     string
  format:   'jplace' | 'newick'
  sha256:   string
  modified: string
  content:  string
  /** Saved view state, as last written by the tree viewer. */
  view:     unknown
}

export const REFERENCE_STEPS = ['align', 'trim', 'tree'] as const
export const PLACEMENT_STEPS = ['align', 'trim', 'place', 'accumulate'] as const
export type PhyloStep = typeof REFERENCE_STEPS[number] | typeof PLACEMENT_STEPS[number]

export interface PlacementRunSpec { run: string; group?: string | null; subgroups?: string[] }

export interface PlacementQueries {
  source: 'taxon' | 'fasta'
  runs: PlacementRunSpec[]
  table: string
  rank: string
  values: string[]
  min_reads: number
}

export interface AlignSettings { strategy?: string; maxiterate?: number; optional_args?: string }
export interface TrimSettings {
  method?: string
  gap_threshold?: number
  conservation?: number | null
  similarity_threshold?: number | null
  residue_overlap?: number | null
  sequence_overlap?: number | null
  optional_args?: string
}

/** The phylogeny section of pipeline.yml, or part of it. */
export interface PhyloSettings {
  threads?: number
  reference?: {
    align?: AlignSettings
    trim?: TrimSettings
    tree?: { model?: string; bootstrap?: string; replicates?: number; optional_args?: string }
  }
  placement?: {
    align?: AlignSettings
    trim?: TrimSettings
    place?: { model?: string; heuristic?: number | null; optional_args?: string }
    accumulate?: { threshold?: number }
  }
}

export type PhyloState = 'new' | 'running' | 'done' | 'failed'

export interface PhyloStepStatus {
  state: 'pending' | 'running' | 'done' | 'current' | 'failed'
  where?: string
  started?: string
  finished?: string
  error?: string
}

/** What a reference tree's and a placement's detail have in common. */
export interface WorkflowDetail {
  status: { state: PhyloState; error?: string; steps?: Partial<Record<PhyloStep, PhyloStepStatus>> } | null
  state: PhyloState
  settings: PhyloSettings
  inherited: PhyloSettings
  /** The server each step would run on, or null for this machine. */
  remote: Partial<Record<PhyloStep, string | null>>
  /** Steps with a QC summary. */
  qc: PhyloStep[]
}

export interface ReferenceTreeDoc {
  id: string
  name: string
  description: string
  overrides: PhyloSettings
  created: string
  modified: string
}

export interface ReferenceTreeSummary {
  id: string
  name: string
  modified: string
  references: number
  built: boolean
  state: PhyloState
}

export interface ReferenceTreeDetail extends WorkflowDetail {
  doc: ReferenceTreeDoc
  summary: ReferenceTreeSummary
  used_by: { study: string; id: string; name: string }[]
}

export interface PlacementDoc {
  id: string
  name: string
  reference: string | null
  queries: PlacementQueries
  overrides: PhyloSettings
  published?: string[]
  created: string
  modified: string
}

export interface PlacementSummary {
  id: string
  name: string
  modified: string
  reference: string | null
  queries: number
  state: PhyloState
}

export interface PlacementDetail extends WorkflowDetail {
  doc: PlacementDoc
  summary: PlacementSummary
}

export interface AlignmentQC {
  kind: 'alignment'
  sequences: number
  columns: number
  gap_fraction: number
  occupancy: number[]
  occupancy_references?: number[]
  occupancy_queries?: number[]
  kept_columns?: number[]
  removed?: string[]
  per_sequence: { name: string; query: boolean; residues: number; kept_residues: number | null; span: [number, number] | null }[]
}

export interface ReadFunnelData {
  stages:  { key: string; label: string }[]
  samples: { sample: string; values: Record<string, number | null> }[]
  /** Steps that should keep every read, and whether they did. */
  checks:  { name: string; ok: boolean; detail: string }[]
}

export type ReportKind = 'figure' | 'table' | 'tree'

export interface ReportItem {
  id:    string
  kind:  ReportKind
  title: string
  file:  string
  added: string
}
