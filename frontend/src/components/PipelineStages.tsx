// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState, useEffect } from 'react'
import type { RunStages, StageStatus, ConfigMap, ConfigSource } from '../api/types'
import { api } from '../api/client'
import { errorMessage } from '../api/errorMessage'
import { useToast } from './Toast'
import { parseUtc, timeAgo } from '../utils/timeago'
import { splitLines } from '../utils/text'
import { SOURCE_COLORS } from './annotationShared'
import styles from './PipelineStages.module.css'

// cdhit config lives under dada2_denoise; merge_taxa is auto-chained and hidden
const VISIBLE_STAGES = ['cutadapt', 'dada2_denoise', 'dada2_classify', 'swarm', 'vsearch'] as const
type VisibleStage = typeof VISIBLE_STAGES[number]

// Config sections include pipeline stages plus non-stage sections like study_design
export type ConfigSection = VisibleStage | 'study_design' | 'analysis' | 'global' | 'remote' | 'phylogeny'

export const STAGE_LABELS: Record<ConfigSection, string> = {
  global:         'Global',
  remote:         'Remote Execution (Bioserver)',
  study_design:   'Study Design',
  cutadapt:       'Primer Trimming',
  dada2_denoise:  'DADA2 Denoising',
  dada2_classify: 'Taxonomy Assignment (DADA2)',
  swarm:          'OTU Clustering (SWARM)',
  vsearch:        'Taxonomy Assignment (VSEARCH)',
  analysis:       'Analysis',
  phylogeny:      'Phylogenetic Placement',
}

export const STAGE_CONFIG_PREFIXES: Record<ConfigSection, string[]> = {
  global:         ['r_threads'],
  remote:         ['remote.'],
  study_design:   ['seed', 'subsample_n', 'pool_children'],
  cutadapt:       ['cutadapt.'],
  dada2_denoise:  ['dada2.file_patterns.', 'dada2.filter_trim.', 'dada2.dada.', 'dada2.merge.', 'dada2.asv.', 'cdhit.'],
  dada2_classify: ['dada2.taxonomy.', 'dada2.output.'],
  swarm:          ['swarm.'],
  vsearch:        ['vsearch.'],
  analysis:       ['analysis.alpha.', 'analysis.beta.', 'analysis.differential.', 'analysis.taxa_bar.'],
  phylogeny:      ['phylogeny.threads', 'phylogeny.reference.align.', 'phylogeny.reference.trim.', 'phylogeny.reference.tree.',
                   'phylogeny.placement.align.', 'phylogeny.placement.trim.', 'phylogeny.placement.place.',
                   'phylogeny.placement.accumulate.'],
}

// Headings for prefixes whose last segment alone would be ambiguous.
export const SECTION_LABELS: Record<string, string> = {
  'phylogeny.reference.align.':      'Reference tree: align (MAFFT)',
  'phylogeny.reference.trim.':       'Reference tree: trim (trimAl)',
  'phylogeny.reference.tree.':       'Reference tree: tree (IQ-TREE)',
  'phylogeny.placement.align.':      'Placement: add queries (MAFFT)',
  'phylogeny.placement.trim.':       'Placement: trim (trimAl)',
  'phylogeny.placement.place.':      'Placement: place (RAxML EPA)',
  'phylogeny.placement.accumulate.': 'Placement: accumulate (gappa)',
}

// Keys the pipeline no longer reads that may linger in a resolved config, such as
// the SWARM pins kept to preserve stage hashes, and deprecated aliases of r_threads
// and remote.*. The editor hides them.
export const HIDDEN_CONFIG_KEYS = new Set<string>([
  'swarm.chimera_check',
  'swarm.min_abundance',
  'swarm.fastq_minovlen',
  'swarm.identity',
  'dada2.filter_trim.max_n',
  'dada2.taxonomy.multithread',
  'dada2.taxonomy.remote.host',
  'dada2.taxonomy.remote.identity_file',
  'dada2.taxonomy.remote.rscript',
  'dada2.taxonomy.remote.staging_dir',
])

const STAGE_ORDER = [...VISIBLE_STAGES]

export const CONFIG_DESCRIPTIONS: Record<string, string> = {
  'r_threads':                      'Worker threads for every DADA2 R stage (learn_errors, denoise, chimera removal, taxonomy). Does not affect stage hashes - changing it never marks a stage stale.',
  'remote.host':                    'SSH destination for offloaded stages, e.g. user@host. Empty runs every stage locally.',
  'remote.rscript':                 'Rscript on the server; a bare name is resolved against the login PATH.',
  'remote.staging_dir':             'Absolute directory on the server for per-run staging. Each stage creates and removes its own subdirectory.',
  'remote.identity_file':           'SSH private key; empty falls back to the agent, then to a password prompt.',
  'remote.threads':                 'Threads to request on the server; empty uses r_threads for DADA2 and phylogeny.threads for the phylogeny steps.',
  'remote.stages':                  'Stages to run on the server. learn_errors and denoise upload the filtered reads; chimera_removal and assign_taxonomy upload only checkpoints. The phylogeny stages upload their FASTA and tree files; trimming and gappa always run here.',
  'remote.tools.mafft':             'MAFFT on the server; a bare name is resolved against the login PATH.',
  'remote.tools.iqtree':            'IQ-TREE on the server, e.g. iqtree, iqtree2 or iqtree3.',
  'remote.tools.raxml':             'RAxML on the server; must be a PTHREADS build.',
  'pool_children':                  'Pool FASTQ files from child directories into a single run. Each child directory becomes a sub-group identified by a prefix.',
  'seed':                           'Global random seed used by stages that require randomness.',
  'subsample_n':                    'Number of samples to randomly select per run for quick config iteration (0 = all). 3 is recommended for testing DADA2 parameters.',
  'cutadapt.primer_pairs':        'Primer pair names to apply. Must match keys in config/primers.yml.',
  'cutadapt.min_length':          'Discard reads shorter than this after trimming (-m).',
  'cutadapt.discard_untrimmed':   'Drop reads where no adapter was found (--discard-untrimmed).',
  'cutadapt.cores':               'Parallel cores; 0 = auto-detect (-j).',
  'cutadapt.quality_cutoff':      "3' quality trimming cutoff; null to disable (-q).",
  'cutadapt.error_rate':          'Max adapter mismatch rate; null = cutadapt default (-e).',
  'cutadapt.overlap':             'Min adapter overlap length; null = cutadapt default (-O).',
  'cutadapt.r1_suffix':           'Filename suffix identifying R1 (forward) reads, e.g. "_R1" or "_1".',
  'cutadapt.r2_suffix':           'Filename suffix identifying R2 (reverse) reads, e.g. "_R2" or "_2".',
  'cutadapt.optional_args':       'Additional flags passed verbatim to cutadapt.',
  'dada2.file_patterns.mode':     'Read mode: paired, forward, or reverse.',
  'dada2.filter_trim.trunc_q':    'Truncate reads at the first base with quality <= this value.',
  'dada2.filter_trim.trunc_len':  '[forward, reverse] truncation lengths; first value used for single-end.',
  'dada2.filter_trim.max_ee':     '[forward, reverse] max expected errors. Reads exceeding this are discarded.',
  'dada2.filter_trim.min_len':    'Discard reads shorter than this after truncation.',
  'dada2.filter_trim.match_ids':  'Require forward/reverse read IDs to match.',
  'dada2.filter_trim.rm_phix':    'Remove reads matching the PhiX genome.',
  'dada2.dada.nbases':            'Number of bases used for error model learning.',
  'dada2.dada.max_consist':       'Max iterations for error model convergence.',
  'dada2.dada.pool_method':       'Sample pooling: none, pseudo (recommended), or true (memory-intensive).',
  'dada2.merge.min_overlap':      'Minimum overlap (bp) required to merge forward/reverse reads.',
  'dada2.merge.max_mismatch':     'Max mismatches allowed in the overlap region.',
  'dada2.merge.trim_overhang':    'Trim overhanging bases beyond the start of the opposite read.',
  'dada2.asv.band_size_min':      'Min ASV length to keep; null to skip length filtering.',
  'dada2.asv.band_size_max':      'Max ASV length to keep.',
  'dada2.asv.denovo_method':      'Chimera removal method: consensus or pooled.',
  'dada2.taxonomy.enabled':       'Set to false to skip DADA2 taxonomy assignment.',
  'dada2.taxonomy.database':      'Reference database key (from config/databases.yml).',
  'dada2.taxonomy.multithread':   'Threads for assignTaxonomy; higher values increase memory use significantly.',
  'dada2.taxonomy.min_boot':      'Minimum bootstrap confidence for taxonomy assignment.',
  'dada2.taxonomy.remote.host':         'SSH user@hostname for remote execution; null = run locally.',
  'dada2.taxonomy.remote.identity_file':'Path to SSH private key; null = password auth.',
  'dada2.taxonomy.remote.rscript':      'Path to Rscript on the remote server.',
  'dada2.taxonomy.remote.staging_dir':  'Absolute path on the remote server for staging files.',
  'dada2.output.seq_table_prefix':'Filename prefix for the sequence table output.',
  'dada2.output.fasta_prefix':    'Filename prefix for the ASV FASTA output.',
  'dada2.output.taxa_prefix':     'Filename prefix for the taxonomy table output.',
  'dada2.verbose':                'Print verbose progress messages during DADA2 execution.',
  'vsearch.enabled':              'Set to false to skip VSEARCH taxonomy assignment entirely.',
  'vsearch.identity':             'Minimum sequence identity threshold (--id).',
  'vsearch.query_cov':            'Minimum fraction of query sequence covered (--query_cov).',
  'vsearch.maxaccepts':           'Stop after this many hits per query; null = vsearch default (--maxaccepts).',
  'vsearch.maxrejects':           'Max rejected candidates per query; null = vsearch default (--maxrejects).',
  'vsearch.strand':               '"plus" or "both"; null = vsearch default (--strand).',
  'vsearch.optional_args':        'Additional flags passed verbatim to vsearch.',
  'cdhit.enabled':                'Run CD-HIT ASV dereplication before OTU clustering.',
  'cdhit.identity':               'Sequence identity threshold (-c); use 1.0 to deduplicate only.',
  'cdhit.threads':                'Worker threads; 0 = all available (-T).',
  'cdhit.optional_args':          'Additional flags passed verbatim to cd-hit-est.',
  'swarm.enabled':                'Set to false to skip OTU clustering entirely.',
  'swarm.differences':            'Max differences between sequences in the same cluster (-d).',
  'swarm.threads':                'Worker threads; 0 = all available (-t).',
  'swarm.optional_args':          'Additional flags passed verbatim to swarm.',
  'analysis.taxa_bar.top_n':      'Collapse taxa below this rank count to "Other".',
  'analysis.taxa_bar.rank':       'Taxonomic rank for bar plots; null = lowest assigned rank.',
  'analysis.taxa_bar.ranks':      'Rank levels to plot; null = auto (last 3 levels).',
  'analysis.taxa_bar.report_ranks':'Ranks to include in the report; null = auto (last 3 levels).',
  'analysis.alpha.normalisation': 'Normalisation before richness, Shannon and Simpson: none (raw counts), rarefy or srs. rarefy and srs leave out samples below the depth.',
  'analysis.alpha.depth':         'Reads per sample for rarefy or srs; 0 uses the smallest non-empty sample.',
  'analysis.alpha.rarefaction_iterations': 'rarefy draws averaged for each alpha diversity value.',
  'analysis.beta.normalisation':  'Preparation before the Bray-Curtis dissimilarity used by NMDS and PERMANOVA: hellinger (square root of relative abundances), none, rarefy or srs.',
  'analysis.beta.depth':          'Reads per sample for rarefy or srs; 0 uses the smallest non-empty sample.',
  'analysis.differential.offset': 'Library-size offset in each negative-binomial model: tss (total reads per sample) or rle (median-of-ratios over taxa present in every sample; refused when there are none).',
  'analysis.differential.min_prevalence': 'Fraction of samples, 0 to 1, in which a taxon must have reads to be tested; taxa below it are listed as filtered.',
  'analysis.alpha.show_points':   'Show individual samples overlaid on alpha comparison boxplots.',
  'analysis.alpha.annotate_significance': 'Annotate alpha comparison panels with Kruskal-Wallis significance.',
  'analysis.alpha.pairwise_brackets': 'Run BH-adjusted pairwise Wilcoxon rank-sum tests between every group pair and draw brackets, including n.s. results.',
  'analysis.alpha.paired_samples': 'Use paired-sample tests by matching sample IDs across groups. Uses paired Wilcoxon for 2 groups and Friedman for 3+ groups.',
  'analysis.alpha.paired_lines':  'Join each sample to itself across groups with a line, matched as in the paired tests. Disables point jitter.',
  'analysis.alpha.significance_test': 'Overall significance test used for alpha comparison annotations.',
  'phylogeny.threads':                        'Threads for MAFFT, IQ-TREE, RAxML and gappa on this machine.',
  'phylogeny.reference.align.strategy':       'MAFFT method for reference trees: auto, localpair (L-INS-i), genafpair (E-INS-i), globalpair (G-INS-i) or 6merpair.',
  'phylogeny.reference.align.maxiterate':     'MAFFT --maxiterate; 0 omits it, auto ignores it.',
  'phylogeny.reference.align.optional_args':  'Additional flags passed verbatim to MAFFT.',
  'phylogeny.reference.trim.method':             'trimAl selection: manual uses -gt, -cons and -st; gappyout, strict, strictplus, automated1, nogaps and noallgaps are trimAl modes.',
  'phylogeny.reference.trim.gap_threshold':      'trimAl -gt: keep columns with residues in at least this fraction of sequences (manual only).',
  'phylogeny.reference.trim.conservation':       'trimAl -cons: keep at least this percentage of columns whatever the thresholds (manual only).',
  'phylogeny.reference.trim.similarity_threshold': 'trimAl -st: minimum average similarity of a kept column (manual only).',
  'phylogeny.reference.trim.residue_overlap':    'trimAl -resoverlap, set with -seqoverlap to remove poorly overlapping sequences.',
  'phylogeny.reference.trim.sequence_overlap':   'trimAl -seqoverlap: percentage of good positions a sequence needs to be kept.',
  'phylogeny.reference.trim.optional_args':      'Additional flags passed verbatim to trimAl.',
  'phylogeny.reference.tree.model':           'IQ-TREE substitution model (-m); MFP picks one with ModelFinder.',
  'phylogeny.reference.tree.bootstrap':       'Branch support: standard (-b) or ultrafast (-bb) bootstrap.',
  'phylogeny.reference.tree.replicates':      'Bootstrap replicates; ultrafast needs at least 1000.',
  'phylogeny.reference.tree.optional_args':   'Additional flags passed verbatim to IQ-TREE.',
  'phylogeny.placement.align.strategy':       'MAFFT method for adding queries with --addfragments.',
  'phylogeny.placement.align.maxiterate':     'MAFFT --maxiterate for adding queries; 0 omits it.',
  'phylogeny.placement.align.optional_args':  'Additional flags passed verbatim to MAFFT.',
  'phylogeny.placement.trim.method':             'trimAl selection: manual uses -gt, -cons and -st; gappyout, strict, strictplus, automated1, nogaps and noallgaps are trimAl modes.',
  'phylogeny.placement.trim.gap_threshold':      'trimAl -gt: keep columns with residues in at least this fraction of sequences (manual only).',
  'phylogeny.placement.trim.conservation':       'trimAl -cons: keep at least this percentage of columns whatever the thresholds (manual only).',
  'phylogeny.placement.trim.similarity_threshold': 'trimAl -st: minimum average similarity of a kept column (manual only).',
  'phylogeny.placement.trim.residue_overlap':    'trimAl -resoverlap, set with -seqoverlap to remove poorly overlapping sequences.',
  'phylogeny.placement.trim.sequence_overlap':   'trimAl -seqoverlap: percentage of good positions a sequence needs to be kept.',
  'phylogeny.placement.trim.optional_args':      'Additional flags passed verbatim to trimAl.',
  'phylogeny.placement.place.model':          'RAxML substitution model (-m) for EPA placement.',
  'phylogeny.placement.place.heuristic':      'EPA heuristic (-G): the fraction of branches tried in full; empty tries all.',
  'phylogeny.placement.place.optional_args':  'Additional flags passed verbatim to RAxML.',
  'phylogeny.placement.accumulate.threshold': 'gappa accumulate: placement mass a branch must gather before a query is assigned to it (0.5 to 1).',
}

export type ConfigType =
  | { kind: 'boolean' }
  | { kind: 'int'; nullable?: boolean }
  | { kind: 'float'; nullable?: boolean; min?: number; max?: number; step?: number }
  | { kind: 'enum'; options: string[] }
  | { kind: 'multiselect'; optionsFrom: string }
  | { kind: 'string_list' }

/** Explicit type hints for config keys that aren't plain strings or arrays. */
export const CONFIG_TYPES: Record<string, ConfigType> = {
  'r_threads':                      { kind: 'int' },
  'remote.threads':                 { kind: 'int', nullable: true },
  'remote.stages':                  { kind: 'multiselect', optionsFrom: 'remote_stages' },
  'pool_children':                  { kind: 'boolean' },
  'seed':                           { kind: 'int' },
  'subsample_n':                    { kind: 'int' },
  'cutadapt.primer_pairs':        { kind: 'multiselect', optionsFrom: 'primers' },
  'cutadapt.min_length':          { kind: 'int' },
  'cutadapt.discard_untrimmed':   { kind: 'boolean' },
  'cutadapt.cores':               { kind: 'int' },
  'cutadapt.quality_cutoff':      { kind: 'int', nullable: true },
  'cutadapt.error_rate':          { kind: 'float', nullable: true, min: 0, max: 1, step: 0.01 },
  'cutadapt.overlap':             { kind: 'int', nullable: true },
  'dada2.file_patterns.mode':     { kind: 'enum', options: ['paired', 'forward', 'reverse'] },
  'dada2.filter_trim.trunc_q':    { kind: 'int' },
  'dada2.filter_trim.min_len':    { kind: 'int' },
  'dada2.filter_trim.match_ids':  { kind: 'boolean' },
  'dada2.filter_trim.rm_phix':    { kind: 'boolean' },
  'dada2.dada.nbases':            { kind: 'int' },
  'dada2.dada.max_consist':       { kind: 'int' },
  'dada2.dada.pool_method':       { kind: 'enum', options: ['none', 'pseudo', 'true'] },
  'dada2.merge.min_overlap':      { kind: 'int' },
  'dada2.merge.max_mismatch':     { kind: 'int' },
  'dada2.merge.trim_overhang':    { kind: 'boolean' },
  'dada2.asv.band_size_min':      { kind: 'int', nullable: true },
  'dada2.asv.band_size_max':      { kind: 'int' },
  'dada2.asv.denovo_method':      { kind: 'enum', options: ['consensus', 'pooled'] },
  'dada2.taxonomy.enabled':       { kind: 'boolean' },
  'dada2.taxonomy.multithread':   { kind: 'int' },
  'dada2.taxonomy.min_boot':      { kind: 'int' },
  'dada2.verbose':                { kind: 'boolean' },
  'vsearch.identity':             { kind: 'float', min: 0, max: 1, step: 0.01 },
  'vsearch.query_cov':            { kind: 'float', min: 0, max: 1, step: 0.01 },
  'vsearch.maxaccepts':           { kind: 'int', nullable: true },
  'vsearch.maxrejects':           { kind: 'int', nullable: true },
  'vsearch.strand':               { kind: 'enum', options: ['plus', 'both'] },
  'vsearch.enabled':              { kind: 'boolean' },
  'cdhit.enabled':                { kind: 'boolean' },
  'cdhit.identity':               { kind: 'float', min: 0, max: 1, step: 0.01 },
  'cdhit.threads':                { kind: 'int' },
  'swarm.enabled':                { kind: 'boolean' },
  'swarm.differences':            { kind: 'int' },
  'swarm.threads':                { kind: 'int' },
  'analysis.taxa_bar.top_n':      { kind: 'int' },
  'analysis.alpha.normalisation': { kind: 'enum', options: ['none', 'rarefy', 'srs'] },
  'analysis.alpha.depth':         { kind: 'int' },
  'analysis.alpha.rarefaction_iterations': { kind: 'int' },
  'analysis.beta.normalisation':  { kind: 'enum', options: ['hellinger', 'none', 'rarefy', 'srs'] },
  'analysis.beta.depth':          { kind: 'int' },
  'analysis.differential.offset': { kind: 'enum', options: ['tss', 'rle'] },
  'analysis.differential.min_prevalence': { kind: 'float', min: 0, max: 1, step: 0.05 },
  'analysis.alpha.show_points':   { kind: 'boolean' },
  'analysis.alpha.annotate_significance': { kind: 'boolean' },
  'analysis.alpha.pairwise_brackets': { kind: 'boolean' },
  'analysis.alpha.paired_samples': { kind: 'boolean' },
  'analysis.alpha.paired_lines':  { kind: 'boolean' },
  'analysis.alpha.significance_test': { kind: 'enum', options: ['kruskal_wallis'] },
  'phylogeny.threads':                        { kind: 'int' },
  'phylogeny.reference.align.strategy':       { kind: 'enum', options: ['auto', 'localpair', 'genafpair', 'globalpair', '6merpair'] },
  'phylogeny.reference.align.maxiterate':     { kind: 'int' },
  'phylogeny.reference.trim.method':          { kind: 'enum', options: ['manual', 'gappyout', 'strict', 'strictplus', 'automated1', 'nogaps', 'noallgaps'] },
  'phylogeny.reference.trim.gap_threshold':   { kind: 'float', min: 0, max: 1, step: 0.01 },
  'phylogeny.reference.trim.conservation':    { kind: 'float', nullable: true, min: 0, max: 100, step: 1 },
  'phylogeny.reference.trim.similarity_threshold': { kind: 'float', nullable: true, min: 0, max: 1, step: 0.001 },
  'phylogeny.reference.trim.residue_overlap': { kind: 'float', nullable: true, min: 0, max: 1, step: 0.01 },
  'phylogeny.reference.trim.sequence_overlap': { kind: 'float', nullable: true, min: 0, max: 100, step: 1 },
  'phylogeny.reference.tree.bootstrap':       { kind: 'enum', options: ['standard', 'ultrafast'] },
  'phylogeny.reference.tree.replicates':      { kind: 'int' },
  'phylogeny.placement.align.strategy':       { kind: 'enum', options: ['auto', 'localpair', 'genafpair', 'globalpair', '6merpair'] },
  'phylogeny.placement.align.maxiterate':     { kind: 'int' },
  'phylogeny.placement.trim.method':          { kind: 'enum', options: ['manual', 'gappyout', 'strict', 'strictplus', 'automated1', 'nogaps', 'noallgaps'] },
  'phylogeny.placement.trim.gap_threshold':   { kind: 'float', min: 0, max: 1, step: 0.01 },
  'phylogeny.placement.trim.conservation':    { kind: 'float', nullable: true, min: 0, max: 100, step: 1 },
  'phylogeny.placement.trim.similarity_threshold': { kind: 'float', nullable: true, min: 0, max: 1, step: 0.001 },
  'phylogeny.placement.trim.residue_overlap': { kind: 'float', nullable: true, min: 0, max: 1, step: 0.01 },
  'phylogeny.placement.trim.sequence_overlap': { kind: 'float', nullable: true, min: 0, max: 100, step: 1 },
  'phylogeny.placement.place.heuristic':      { kind: 'float', nullable: true, min: 0, max: 1, step: 0.01 },
  'phylogeny.placement.accumulate.threshold': { kind: 'float', min: 0.5, max: 1, step: 0.01 },
}


interface Props {
  stages:   RunStages
  onRun?:   (stage: string) => void | Promise<unknown>
  disabled?: boolean
  configMap?: ConfigMap | null
  study?:    string
  run?:      string
  group?:    string
  onConfigChanged?: () => void
}

function StatusDot({ status }: { status: StageStatus }) {
  const label = status.replace('_', ' ')
  return <span className={`${styles.dot} ${styles[status]}`} title={label} role="img" aria-label={label} />
}

export function PipelineStages({ stages, onRun, disabled, configMap, study, run, group, onConfigChanged }: Props) {
  const [expanded, setExpanded] = useState<string | null>(null)
  const [pending, setPending]   = useState<Set<string>>(new Set())

  useEffect(() => {
    setPending(prev => {
      if (prev.size === 0) return prev
      const next = new Set(prev)
      for (const key of prev) {
        const status = stages?.[key as keyof RunStages]?.status
        if (status === 'running' || status === 'complete' || status === 'stale' || status === 'not_started')
          next.delete(key)
      }
      return next.size === prev.size ? prev : next
    })
  }, [stages])

  return (
    <div className={styles.grid}>
      {STAGE_ORDER.map(key => {
        const info = stages?.[key]
        const status   = info?.status ?? 'not_started'
        const last_run = info?.last_run ?? null
        const isExpanded = expanded === key
        const hasConfig = configMap && STAGE_CONFIG_PREFIXES[key].some(
          prefix => Object.keys(configMap).some(k => k.startsWith(prefix))
        )
        // Count how many config keys are overridden at run level for this stage
        const runOverrides = configMap ? STAGE_CONFIG_PREFIXES[key].reduce((n, prefix) =>
          n + Object.entries(configMap).filter(([k, { source }]) => k.startsWith(prefix) && source === 'run').length, 0) : 0
        return (
          <div key={key}>
            <div className={`${styles.row} ${styles[status]}`}>
              <StatusDot status={status} />
              <button
                type="button"
                className={`toggle-btn ${styles.label}`}
                style={{ cursor: hasConfig ? 'pointer' : 'default', display: 'inline', width: 'auto' }}
                aria-expanded={hasConfig ? isExpanded : undefined}
                disabled={!hasConfig}
                onClick={() => hasConfig && setExpanded(isExpanded ? null : key)}
              >
                {hasConfig && <span aria-hidden="true" style={{ fontSize: '.8rem', marginRight: 6, opacity: .65 }}>{isExpanded ? '▾' : '▸'}</span>}
                {STAGE_LABELS[key]}
                {runOverrides > 0 && (
                  <span style={{ marginLeft: 6, fontSize: '.68rem', fontWeight: 600, color: 'var(--color-primary)', verticalAlign: 'middle' }}
                    title={`${runOverrides} run-level override${runOverrides > 1 ? 's' : ''}`}>
                    {runOverrides} override{runOverrides > 1 ? 's' : ''}
                  </span>
                )}
              </button>
              <span className={styles.ts} title={last_run ? parseUtc(last_run).toLocaleString() : ''}>{last_run ? timeAgo(last_run) : '-'}</span>
              {onRun && status !== 'disabled' && (
                <button
                  className={styles.run}
                  disabled={disabled || status === 'running' || pending.has(key)}
                  onClick={() => {
                    setPending(prev => new Set(prev).add(key))
                    // A rejected request never changes `stages`, so clear the pending mark here.
                    Promise.resolve(onRun(key)).catch(() => setPending(prev => {
                      const next = new Set(prev)
                      next.delete(key)
                      return next
                    }))
                  }}
                >
                  {status === 'running' || pending.has(key) ? '…' : 'Run'}
                </button>
              )}
            </div>
            {isExpanded && configMap && study && run && onConfigChanged && (
              <StageConfig
                configMap={configMap}
                prefixes={STAGE_CONFIG_PREFIXES[key]}
                study={study}
                run={run}
                group={group}
                onConfigChanged={onConfigChanged}
              />
            )}
          </div>
        )
      })}
    </div>
  )
}

export function StageConfig({ configMap, prefixes, labels, study, run, group, onConfigChanged, patchFn, deleteFn, sourceLevel, overrides }: {
  configMap: ConfigMap
  prefixes: string[]
  /** Headings by prefix, in place of SECTION_LABELS. */
  labels?: Record<string, string>
  study: string
  run: string
  group?: string
  onConfigChanged: () => void
  patchFn?: (study: string, run: string, body: Record<string, unknown>, group?: string) => Promise<ConfigMap>
  deleteFn?: (study: string, run: string, key: string, group?: string) => Promise<ConfigMap>
  sourceLevel?: ConfigSource
  overrides?: Record<string, string[]> | null
}) {
  // Group entries by section prefix for nice headers
  const sections: { label: string | null; entries: { dottedKey: string; leafKey: string; value: unknown; source: ConfigSource }[] }[] = []

  for (const prefix of prefixes) {
    const isSectionPrefix = prefix.endsWith('.')
    const sectionEntries = Object.entries(configMap)
      .filter(([k]) => !HIDDEN_CONFIG_KEYS.has(k) && (isSectionPrefix ? k.startsWith(prefix) : k === prefix))
      .map(([k, { value, source }]) => ({
        dottedKey: k,
        leafKey: isSectionPrefix ? k.slice(prefix.length) : k,
        value,
        source,
      }))
      .sort((a, b) => a.leafKey.localeCompare(b.leafKey))
    if (sectionEntries.length > 0) {
      const label = !isSectionPrefix ? null
        : labels?.[prefix] ?? SECTION_LABELS[prefix] ?? (prefix.replace(/\.$/, '').split('.').pop() ?? prefix)
            .replace(/_/g, ' ')
            .replace(/\b\w/g, c => c.toUpperCase())
      sections.push({ label, entries: sectionEntries })
    }
  }

  if (sections.length === 0) return null

  return (
    <div className={styles.configPanel}>
      {sections.map(section => (
        <div key={section.label ?? section.entries.map(e => e.dottedKey).join('|')}>
          {section.label && (
            <div style={{ fontWeight: 600, fontSize: '.78rem', marginBottom: 2, marginTop: 4 }}>{section.label}</div>
          )}
          {section.entries.map(e => (
            <StageConfigField
              key={e.dottedKey}
              dottedKey={e.dottedKey}
              leafKey={e.leafKey}
              value={e.value}
              source={e.source}
              study={study}
              run={run}
              group={group}
              onChanged={onConfigChanged}
              patchFn={patchFn}
              deleteFn={deleteFn}
              sourceLevel={sourceLevel}
              overrides={overrides?.[e.dottedKey]}
            />
          ))}
        </div>
      ))}
    </div>
  )
}

function StageConfigField({ dottedKey, leafKey, value, source, study, run, group, onChanged, patchFn, deleteFn, sourceLevel = 'run', overrides }: {
  dottedKey: string
  leafKey: string
  value: unknown
  source: ConfigSource
  study: string
  run: string
  group?: string
  onChanged: () => void
  patchFn?: (study: string, run: string, body: Record<string, unknown>, group?: string) => Promise<ConfigMap>
  deleteFn?: (study: string, run: string, key: string, group?: string) => Promise<ConfigMap>
  sourceLevel?: ConfigSource
  overrides?: string[]
}) {
  const [editing, setEditing] = useState(false)
  const [draft, setDraft] = useState('')
  const [saving, setSaving] = useState(false)
  const toast = useToast()

  const tooltip = CONFIG_DESCRIPTIONS[dottedKey]
  const typeHint = CONFIG_TYPES[dottedKey]

  const displayValue = Array.isArray(value) ? JSON.stringify(value)
    : value === null || value === undefined ? 'null'
    : typeof value === 'boolean' ? (value ? 'true' : 'false')
    : typeof value === 'object' ? JSON.stringify(value)
    : String(value)

  const startEdit = () => {
    // Booleans, enums, and multiselects use inline controls, no draft needed
    if (typeHint?.kind === 'boolean' || typeHint?.kind === 'enum' || typeHint?.kind === 'multiselect') return
    if (typeHint?.kind === 'string_list') {
      const items = Array.isArray(value) ? (value as string[]) : []
      setDraft(items.join('\n'))
    } else {
      setDraft(typeof value === 'string' ? value : JSON.stringify(value))
    }
    setEditing(true)
  }

  const saveValue = async (newValue: unknown) => {
    setSaving(true)
    try {
      const patch = patchFn ?? api.config.patchRun
      await patch(study, run, { [dottedKey]: newValue }, group)
      setEditing(false)
      onChanged()
    } catch (err) {
      toast.error('Failed to save: ' + (errorMessage(err)))
    } finally {
      setSaving(false)
    }
  }

  const save = async () => {
    let parsed: unknown
    if (typeHint?.kind === 'string_list') {
      parsed = splitLines(draft)
    } else {
      try { parsed = JSON.parse(draft) } catch { parsed = draft }
    }
    await saveValue(parsed)
  }

  const remove = async () => {
    if (!window.confirm(`Remove the ${sourceLevel} override of ${dottedKey}?`)) return
    setSaving(true)
    try {
      const del = deleteFn ?? api.config.deleteRun
      await del(study, run, dottedKey, group)
      setEditing(false)
      onChanged()
    } catch (err) {
      toast.error('Failed to remove: ' + (errorMessage(err)))
    } finally {
      setSaving(false)
    }
  }

  const keyLabel = (leafKey || dottedKey).replace(/_/g, ' ')
  // The rows open an editor on click, so Enter and Space do the same.
  const activate = (e: React.KeyboardEvent, fn: () => void) => {
    if (e.target !== e.currentTarget) return
    if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); fn() }
  }
  const labelEl = (
    <div style={{ minWidth: 120, fontSize: '.78rem' }}>
      {keyLabel}
      {overrides && overrides.length > 0 && (
        <div style={{ fontSize: '.68rem', color: 'var(--color-warn)', lineHeight: 1.2 }}
          title={overrides.join(', ')}>
          {overrides.length} override{overrides.length > 1 ? 's' : ''}
        </div>
      )}
    </div>
  )
  const sourceEl = <span style={{ fontSize: '.68rem', fontWeight: 600, color: SOURCE_COLORS[source], textTransform: 'uppercase', whiteSpace: 'nowrap' }}>{source}</span>
  const removeBtn = source === sourceLevel && (
    <button className="btn" style={{ padding: '0 4px', fontSize: '.68rem', lineHeight: 1 }} title={`Remove ${sourceLevel} override`}
      aria-label={`Remove ${sourceLevel} override of ${keyLabel}`}
      onClick={e => { e.stopPropagation(); remove() }}>&times;</button>
  )

  if (typeHint?.kind === 'boolean') {
    return (
      <div style={{ display: 'flex', gap: 8, alignItems: 'center', padding: '2px 0 2px 8px' }} title={tooltip}>
        {labelEl}
        <input type="checkbox" checked={!!value} disabled={saving}
          onChange={e => saveValue(e.target.checked)}
          style={{ accentColor: 'var(--color-primary)' }} />
        <span style={{ fontFamily: 'monospace', fontSize: '.78rem', flex: 1 }}>{value ? 'true' : 'false'}</span>
        {sourceEl}{removeBtn}
      </div>
    )
  }

  if (typeHint?.kind === 'enum') {
    const nullable = value === null || value === undefined
    return (
      <div style={{ display: 'flex', gap: 8, alignItems: 'center', padding: '2px 0 2px 8px' }} title={tooltip}>
        {labelEl}
        <select
          value={nullable ? '' : String(value)}
          disabled={saving}
          onChange={e => saveValue(e.target.value || null)}
          style={{ fontFamily: 'monospace', fontSize: '.78rem', padding: '1px 4px', border: '1px solid var(--color-border)', borderRadius: 3, background: 'var(--color-bg)' }}
        >
          {nullable && <option value="">null</option>}
          {typeHint.options.map(o => <option key={o} value={o}>{o}</option>)}
        </select>
        <div style={{ flex: 1 }} />
        {sourceEl}{removeBtn}
      </div>
    )
  }

  if (typeHint?.kind === 'multiselect') {
    return <MultiSelectField
      value={value} tooltip={tooltip} optionsFrom={typeHint.optionsFrom}
      saving={saving} saveValue={saveValue}
      labelEl={labelEl} sourceEl={sourceEl} removeBtn={removeBtn} keyLabel={keyLabel}
    />
  }

  if (typeHint?.kind === 'string_list') {
    const items = Array.isArray(value) ? (value as string[]) : []
    if (editing) {
      return (
        <div style={{ display: 'flex', gap: 6, alignItems: 'flex-start', padding: '2px 0 2px 8px' }} title={tooltip}>
          {labelEl}
          <textarea
            aria-label={keyLabel}
            style={{ flex: 1, fontFamily: 'monospace', fontSize: '.78rem', padding: '2px 6px', border: '1px solid var(--color-border)', borderRadius: 3, resize: 'vertical', minHeight: 60, background: 'var(--color-bg)', color: 'inherit' }}
            value={draft}
            onChange={e => setDraft(e.target.value)}
            autoFocus disabled={saving}
            placeholder="one entry per line"
          />
          <div style={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
            <button className="btn" style={{ padding: '1px 6px', fontSize: '.72rem' }} onClick={save} disabled={saving}>Save</button>
            <button className="btn" style={{ padding: '1px 6px', fontSize: '.72rem' }} onClick={() => setEditing(false)} disabled={saving}>Cancel</button>
          </div>
        </div>
      )
    }
    const preview = items.length === 0 ? '(empty)'
      : items.length <= 3 ? items.join(', ')
      : `${items.slice(0, 3).join(', ')} ... (${items.length})`
    return (
      <div
        role="button" tabIndex={0} aria-label={`Edit ${keyLabel}`}
        style={{ display: 'flex', gap: 8, alignItems: 'center', padding: '2px 0 2px 8px', cursor: 'pointer' }}
        onClick={startEdit}
        onKeyDown={e => activate(e, startEdit)}
        title={tooltip ?? 'Click to edit'}
      >
        {labelEl}
        <div style={{ fontFamily: 'monospace', fontSize: '.78rem', flex: 1, color: items.length === 0 ? 'var(--color-muted-fg)' : undefined }}>{preview}</div>
        {sourceEl}{removeBtn}
      </div>
    )
  }

  if (editing) {
    const isNumeric = typeHint?.kind === 'int' || typeHint?.kind === 'float'
    const nullable = isNumeric && typeHint.nullable
    return (
      <div style={{ display: 'flex', gap: 6, alignItems: 'center', padding: '2px 0 2px 8px' }} title={tooltip}>
        {labelEl}
        {isNumeric ? (
          <input
            type="number"
            aria-label={keyLabel}
            style={{ width: 100, fontFamily: 'monospace', fontSize: '.78rem', padding: '2px 6px', border: '1px solid var(--color-border)', borderRadius: 3 }}
            value={draft}
            step={typeHint.kind === 'float' ? (typeHint.step ?? 0.01) : 1}
            min={typeHint.kind === 'float' ? typeHint.min : undefined}
            max={typeHint.kind === 'float' ? typeHint.max : undefined}
            onChange={e => setDraft(e.target.value)}
            onKeyDown={e => { if (e.key === 'Enter') save(); if (e.key === 'Escape') setEditing(false) }}
            autoFocus disabled={saving}
          />
        ) : (
          <input
            aria-label={keyLabel}
            style={{ flex: 1, fontFamily: 'monospace', fontSize: '.78rem', padding: '2px 6px', border: '1px solid var(--color-border)', borderRadius: 3 }}
            value={draft}
            onChange={e => setDraft(e.target.value)}
            onKeyDown={e => { if (e.key === 'Enter') save(); if (e.key === 'Escape') setEditing(false) }}
            autoFocus disabled={saving}
          />
        )}
        {nullable && (
          <button className="btn" style={{ padding: '1px 6px', fontSize: '.72rem' }}
            onClick={() => saveValue(null)} disabled={saving}>null</button>
        )}
        <button className="btn" style={{ padding: '1px 6px', fontSize: '.72rem' }} onClick={save} disabled={saving}>Save</button>
        <button className="btn" style={{ padding: '1px 6px', fontSize: '.72rem' }} onClick={() => setEditing(false)} disabled={saving}>Cancel</button>
      </div>
    )
  }

  return (
    <div
      role="button" tabIndex={0} aria-label={`Edit ${keyLabel}`}
      style={{ display: 'flex', gap: 8, alignItems: 'center', padding: '2px 0 2px 8px', cursor: 'pointer' }}
      onClick={startEdit}
      onKeyDown={e => activate(e, startEdit)}
      title={tooltip ?? 'Click to edit'}
    >
      {labelEl}
      <div style={{ fontFamily: 'monospace', fontSize: '.78rem', flex: 1 }}>{displayValue}</div>
      {sourceEl}{removeBtn}
    </div>
  )
}

const MULTISELECT_FETCHERS: Record<string, () => Promise<string[]>> = {
  primers: () => api.primers.list(),
  // The offloadable stages are fixed by the pipeline. This list mirrors
  // Validation.REMOTE_STAGES, which refuses anything else at the write gate.
  remote_stages: () => Promise.resolve(
    ['learn_errors', 'denoise', 'chimera_removal', 'assign_taxonomy',
     'phylogeny_align', 'phylogeny_tree', 'phylogeny_add', 'phylogeny_place']),
}

function MultiSelectField({ value, tooltip, optionsFrom, saving, saveValue,
  labelEl, sourceEl, removeBtn, keyLabel,
}: {
  value: unknown
  tooltip?: string; optionsFrom: string; saving: boolean
  saveValue: (v: unknown) => Promise<void>
  labelEl: React.ReactNode; sourceEl: React.ReactNode; removeBtn: React.ReactNode
  keyLabel: string
}) {
  const [options, setOptions] = useState<string[] | null>(null)
  const [open, setOpen] = useState(false)
  const selected = Array.isArray(value) ? value.map(String) : []

  useEffect(() => {
    const fetcher = MULTISELECT_FETCHERS[optionsFrom]
    if (fetcher) fetcher().then(setOptions).catch(() => setOptions([]))
  }, [optionsFrom])

  const toggle = (opt: string) => {
    const next = selected.includes(opt)
      ? selected.filter(s => s !== opt)
      : [...selected, opt]
    saveValue(next)
  }

  return (
    <div style={{ padding: '2px 0 2px 8px' }} title={tooltip}>
      <div style={{ display: 'flex', gap: 8, alignItems: 'center' }}>
        {labelEl}
        <div
          role="button" tabIndex={0} aria-expanded={open} aria-label={`Choose ${keyLabel}`}
          style={{ fontFamily: 'monospace', fontSize: '.78rem', flex: 1, cursor: 'pointer', padding: '1px 4px', border: '1px solid var(--color-border)', borderRadius: 3, minHeight: 22, display: 'flex', flexWrap: 'wrap', gap: 4, alignItems: 'center' }}
          onClick={() => setOpen(!open)}
          onKeyDown={e => {
            if (e.target === e.currentTarget && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); setOpen(!open) }
          }}
        >
          {selected.length > 0
            ? selected.map(s => (
                <span key={s} style={{ background: 'var(--color-primary)', color: '#fff', borderRadius: 3, padding: '0 5px', fontSize: '.72rem', lineHeight: '18px' }}>
                  {s}
                </span>
              ))
            : <span style={{ color: 'var(--color-muted-fg)' }}>none</span>}
          <span style={{ marginLeft: 'auto', fontSize: '.68rem', opacity: .5 }}>{open ? 'Hide' : 'Show'}</span>
        </div>
        {sourceEl}{removeBtn}
      </div>
      {open && options && (
        <div style={{ marginLeft: 128, marginTop: 4, border: '1px solid var(--color-border)', borderRadius: 4, padding: 6, maxWidth: 300 }}>
          {options.map(opt => (
            <label key={opt} style={{ display: 'flex', gap: 6, alignItems: 'center', fontSize: '.78rem', cursor: 'pointer', padding: '2px 0' }}>
              <input type="checkbox" checked={selected.includes(opt)} disabled={saving}
                onChange={() => toggle(opt)} style={{ accentColor: 'var(--color-primary)' }} />
              {opt}
            </label>
          ))}
          {options.length === 0 && <div style={{ fontSize: '.75rem', color: 'var(--color-muted-fg)' }}>No options available</div>}
        </div>
      )}
    </div>
  )
}
