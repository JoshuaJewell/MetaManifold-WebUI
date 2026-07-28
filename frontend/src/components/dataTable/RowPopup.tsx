// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { forwardRef } from 'react'
import type { RowPopupData } from '../DataTable'
import styles from '../DataTable.module.css'
import { blastUrl, flashCopy } from './copy'

interface Props {
  data:         RowPopupData | null
  loading:      boolean
  hiddenCols:   Set<string>
  pos:          { x: number; y: number }
  onMouseEnter: () => void
  onMouseLeave: () => void
}

// The ASV members of a hovered row.
export const RowPopup = forwardRef<HTMLDivElement, Props>(function RowPopup(
  { data, loading, hiddenCols, pos, onMouseEnter, onMouseLeave }, ref,
) {
  return (
    <div
      ref={ref}
      className={styles.popup}
      style={{ left: pos.x, top: pos.y }}
      onMouseEnter={onMouseEnter}
      onMouseLeave={onMouseLeave}
    >
      {loading && <div className={styles.popupTitle}>Loading…</div>}
      {!loading && data && data.rows.length > 0 && (() => {
        const popupCols = data.columns.filter(c => !hiddenCols.has(c))
        const popupHasSeq = popupCols.includes('sequence')
        return (
          <>
            <div className={styles.popupTitle}>
              ASV members ({data.rows.length})
            </div>
            <div className={styles.popupScroll}>
              <table className={styles.popupTable}>
                <thead>
                  <tr>
                    {popupCols.map(c => (
                      <th key={c}>{c}</th>
                    ))}
                    {popupHasSeq && <th style={{ width: 50 }}></th>}
                  </tr>
                </thead>
                <tbody>
                  {data.rows.map((r, i) => (
                    <tr key={i}>
                      {popupCols.map(c => {
                        const text = String(r[c] ?? '')
                        return (
                          <td key={c}
                            onClick={flashCopy(text)}
                            title="Click to copy"
                          >{text}</td>
                        )
                      })}
                      {popupHasSeq && (
                        <td className={styles.blastCell}>
                          <a
                            href={blastUrl(String(r['sequence'] ?? ''))}
                            target="_blank"
                            rel="noopener noreferrer"
                            className={styles.blastLink}
                            title="Search this sequence on NCBI BLAST"
                          >BLAST</a>
                        </td>
                      )}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </>
        )
      })()}
    </div>
  )
})
