// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).

/**
 * Turn the contents of a multi-line textarea into a list of values: one per
 * line, trimmed, with blank lines dropped. The shape every "paste a list
 * freely" field in the app expects.
 */
export function splitLines(text: string): string[] {
  return text.split('\n').map(v => v.trim()).filter(v => v.length > 0)
}

/**
 * Strip the pseudo-HTML markup Plotly accepts in titles and axis labels (<b>,
 * <br>, <sub> and friends), leaving the bare text for measuring or for a
 * filename. A single pass is not enough: replacing "<scr<script>ipt>" once
 * removes the inner tag and leaves the halves joined into another, so the
 * strip repeats until the text stops changing.
 */
export function stripHtmlTags(text: string): string {
  let stripped = text
  let previous = ''
  do {
    previous = stripped
    stripped = stripped.replace(/<[^>]+>/g, '')
  } while (stripped !== previous)
  return stripped
}
