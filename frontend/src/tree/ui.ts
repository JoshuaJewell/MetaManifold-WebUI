// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).


export const PLACE_SELECTED = '#c2255c'
export const SELECT_COLOUR = '#228be6'

export const fmt = (x: number, digits = 4) => Number.isFinite(x) ? Number(x.toPrecision(digits)).toString() : ''

export const fmtLwr = (x: number | null) => x === null ? 'Missing' : fmt(x, 3)

// Lower-case words that follow a genus without being part of its name.
const NOT_EPITHET = new Set(['isolate', 'isolates', 'genotype', 'clade', 'group', 'and', 'or', 'incl', 'sensu', 'lato', 'stricto'])

/**
 * A label split into runs, with genus and species names marked italic: a capitalised word
 * or abbreviated genus ("G.") and up to two lower-case epithets after it. Accessions,
 * subtypes such as ST5 or AII, "sp." and words like "isolates" stay upright, as do the
 * capitalised words in `exceptions` (phylogeny.viewer.italics.exceptions), such as "Outgroup".
 */
export function taxonRuns(label: string, exceptions: ReadonlySet<string> = new Set()): { text: string; italic: boolean }[] {
  const runs: { text: string; italic: boolean }[] = []
  let epithets = 0
  for (const part of label.split(/(\s+)/)) {
    let italic: boolean
    if (!part || /^\s+$/.test(part)) italic = runs.length > 0 && runs[runs.length - 1].italic
    else if (/^[A-Z][a-z]+$|^[A-Z][a-z]?\.$/.test(part) && !exceptions.has(part)) { italic = true; epithets = 2 }
    else if (epithets > 0 && /^[a-z][a-z-]*$/.test(part) && !NOT_EPITHET.has(part)) { italic = true; epithets-- }
    else { italic = false; epithets = 0 }
    const last = runs[runs.length - 1]
    if (last && last.italic === italic) last.text += part
    else runs.push({ text: part, italic })
  }
  // A space between an italic name and upright text belongs to the upright run.
  for (let k = 0; k < runs.length - 1; k++) {
    const m = runs[k].italic && !runs[k + 1].italic ? runs[k].text.match(/\s+$/) : null
    if (m) { runs[k].text = runs[k].text.slice(0, -m[0].length); runs[k + 1].text = m[0] + runs[k + 1].text }
  }
  return runs.filter(r => r.text)
}
