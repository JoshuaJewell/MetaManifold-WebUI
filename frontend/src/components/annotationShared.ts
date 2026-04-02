// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).

import type { CSSProperties } from 'react'
import type { AnnotationSource } from '../api/types'

export type ContamStatus = 'unassigned' | 'yes' | 'no'

export const SOURCES: AnnotationSource[] = ['VSEARCH', 'DADA2']

export const RANK_COL: Record<string, Record<AnnotationSource, string>> = {
  species:    { VSEARCH: 'Species',    DADA2: 'Species_dada2' },
  genus:      { VSEARCH: 'Genus',      DADA2: 'Genus_dada2' },
  family:     { VSEARCH: 'Family',     DADA2: 'Family_dada2' },
  order:      { VSEARCH: 'Order',      DADA2: 'Order_dada2' },
  class:      { VSEARCH: 'Class',      DADA2: 'Class_dada2' },
  division:   { VSEARCH: 'Division',   DADA2: 'Division_dada2' },
  supergroup: { VSEARCH: 'Supergroup', DADA2: 'Supergroup_dada2' },
}

export const CONTAM_STYLE: Record<ContamStatus, CSSProperties> = {
  unassigned: { color: 'var(--color-muted-fg)' },
  yes:        { color: '#c92a2a', fontWeight: 600 },
  no:         { color: '#2b8a3e', fontWeight: 600 },
}

export const CONTAM_FILTER_FIELDS: Array<{ key: string; label: string }> = [
  { key: 'function', label: 'Function' },
  { key: 'detailed_function', label: 'Detailed function' },
  { key: 'associated_organism', label: 'Assoc. organism' },
  { key: 'associated_material', label: 'Assoc. material' },
  { key: 'environment', label: 'Environment' },
]

export const FUNCDB_FIELDS = [
  { key: 'Domain',                   label: 'Domain',                   required: false },
  { key: 'supergroup',               label: 'Supergroup',               required: false },
  { key: 'division',                 label: 'Division',                 required: false },
  { key: 'class',                    label: 'Class',                    required: false },
  { key: 'order',                    label: 'Order',                    required: false },
  { key: 'family',                   label: 'Family',                   required: false },
  { key: 'Genus',                    label: 'Genus',                    required: false },
  { key: 'Species',                  label: 'Species',                  required: false },
  { key: 'Assignment_level',         label: 'Assignment level',         required: false },
  { key: 'Function',                 label: 'Function',                 required: true },
  { key: 'Detailed_function',        label: 'Detailed function',        required: false },
  { key: 'Associated_organism',      label: 'Associated organism',      required: false },
  { key: 'Associated_material',      label: 'Associated material',      required: false },
  { key: 'Environment',              label: 'Environment',              required: false },
  { key: 'Potential_human_pathogen', label: 'Potential human pathogen', required: false },
  { key: 'Comment',                  label: 'Comment',                  required: false },
  { key: 'Reference',                label: 'Reference',                required: false },
]

export function prefillFromRow(row: Record<string, unknown>, source: AnnotationSource): Record<string, string> {
  const str = (v: unknown) => (v == null ? '' : String(v))
  const result: Record<string, string> = {}

  const taxMap: Array<[string, string]> = source === 'VSEARCH'
    ? [['Domain', 'Domain'], ['Supergroup', 'supergroup'], ['Division', 'division'], ['Class', 'class'],
       ['Order', 'order'], ['Family', 'family'], ['Genus', 'Genus'], ['Species', 'Species']]
    : [['Supergroup_dada2', 'supergroup'], ['Division_dada2', 'division'], ['Class_dada2', 'class'],
       ['Order_dada2', 'order'], ['Family_dada2', 'family'], ['Genus_dada2', 'Genus'], ['Species_dada2', 'Species']]

  for (const [src, dest] of taxMap) {
    const v = str(row[src])
    if (v) result[dest] = v
  }

  const valueMap: Array<[string, string]> = [
    ['funcdb_function', 'Function'],
    ['funcdb_detailed_function', 'Detailed_function'],
    ['funcdb_associated_organism', 'Associated_organism'],
    ['funcdb_associated_material', 'Associated_material'],
    ['funcdb_environment', 'Environment'],
    ['funcdb_potential_human_pathogen', 'Potential_human_pathogen'],
    ['funcdb_comment', 'Comment'],
    ['funcdb_reference', 'Reference'],
  ]

  for (const [src, dest] of valueMap) {
    const v = str(row[src])
    if (v) result[dest] = v
  }

  const rank = str(row.funcdb_match_rank)
  if (rank && rank !== 'unmatched') result.Assignment_level = rank

  return result
}
