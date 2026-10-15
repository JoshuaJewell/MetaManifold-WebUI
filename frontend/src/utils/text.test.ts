// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { describe, expect, test } from 'bun:test'
import { splitLines, stripHtmlTags } from './text'

describe('stripHtmlTags', () => {
  test('removes the markup Plotly allows in a title', () => {
    expect(stripHtmlTags('CO<sub>2</sub> by <b>run</b>')).toBe('CO2 by run')
  })

  test('a tag hiding inside another does not survive', () => {
    expect(stripHtmlTags('<scr<script>ipt>alert(1)</script>')).not.toMatch(/<[^>]+>/)
  })

  test('text without markup, including a bare >, is unchanged', () => {
    expect(stripHtmlTags('p > 0.05 for all')).toBe('p > 0.05 for all')
  })
})

describe('splitLines', () => {
  test('one trimmed value per line, blanks dropped', () => {
    expect(splitLines(' a\n\n b \nc')).toEqual(['a', 'b', 'c'])
  })
})
