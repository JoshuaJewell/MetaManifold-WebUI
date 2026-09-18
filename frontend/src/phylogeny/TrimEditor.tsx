// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useState } from 'react'
import { errorMessage } from '../api/errorMessage'
import type { AlignmentQC, PhyloSettings, TrimSettings } from '../api/types'
import { AlignmentQCView } from './qc'
import styles from './Phylo.module.css'

/** Runs trimAl with the current trim settings on the current alignment, without saving anything. */
export function TrimPreviewButton({ workflow, overrides, inherited, preview, aligned, onPreview }: {
  workflow: 'reference' | 'placement'
  overrides: PhyloSettings
  inherited: PhyloSettings
  preview: (trim: TrimSettings) => Promise<AlignmentQC>
  /** Whether the align step has run, so there is something to trim. */
  aligned: boolean
  onPreview: (p: TrimPreview) => void
}) {
  const [error, setError] = useState<string | null>(null)
  const [busy, setBusy] = useState(false)

  const run = async () => {
    const trim: TrimSettings = {
      ...(inherited[workflow]?.trim ?? {}),
      ...(overrides[workflow]?.trim ?? {}),
    }
    setBusy(true); setError(null)
    try {
      onPreview({ qc: await preview(trim), gt: trim.method === 'manual' ? trim.gap_threshold ?? null : null })
    } catch (e) { setError(errorMessage(e)) }
    setBusy(false)
  }

  return (
    <>
      <div className={styles.row}>
        <button className="btn btn-sm" disabled={busy || !aligned} onClick={run}
          title={aligned ? undefined : 'Run the align step first'}>{busy ? 'Trimming…' : 'Preview trimming'}</button>
      </div>
      {error && <p className="error-msg">{error}</p>}
    </>
  )
}

export interface TrimPreview { qc: AlignmentQC; gt: number | null }

/** A trimming preview, across the page. */
export function TrimPreviewPanel({ preview, loadAlignment, onClose }: {
  preview: TrimPreview
  loadAlignment: () => Promise<string | null>
  onClose: () => void
}) {
  return (
    <div className={`${styles.section} ${styles.wide}`}>
      <div className={styles.heading}>
        Trimming preview
        <span className={styles.spacer} />
        <button className={styles.linkBtn} onClick={onClose}>Close</button>
      </div>
      <div className={styles.muted}>These settings applied to the current alignment. Run the workflow to use them.</div>
      <AlignmentQCView qc={preview.qc} threshold={preview.gt} loadAlignment={loadAlignment} />
    </div>
  )
}
