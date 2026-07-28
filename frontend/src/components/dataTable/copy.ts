// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import styles from '../DataTable.module.css'

export const blastUrl = (seq: string) =>
  `https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&DATABASE=nt&CMD=Put&ENTREZ_QUERY=NOT+uncultured+organism%5Borganism%5D+NOT+environmental+sample%5Borganism%5D&QUERY=${encodeURIComponent(seq)}`

// Copy text to the clipboard and briefly flash the clicked element as feedback.
// The clipboard API is missing outside secure contexts (plain HTTP to a non-localhost host).
export const flashCopy = (text: string) => (e: React.MouseEvent<HTMLElement>) => {
  if (!navigator.clipboard) return
  navigator.clipboard.writeText(text).catch(() => {})
  const el = e.currentTarget
  el.classList.remove(styles.copied)
  void el.offsetWidth
  el.classList.add(styles.copied)
}
