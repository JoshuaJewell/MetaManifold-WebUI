// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
// SPDX-License-Identifier: MPL-2.0

import { expect, test } from 'bun:test'
import { mkdtempSync, writeFileSync, unlinkSync, rmdirSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { csvField } from './differentialTable'

// Independent character-state reader: doubled quotes are decoded only inside
// a quoted field, and comment markers are active only outside quotes.
function readRecord(record: string): string[] {
  const fields: string[] = []
  let value = '', quoted = false
  for (let i = 0; i < record.length; i++) {
    const c = record[i]
    if (c === '"') {
      if (quoted && record[i + 1] === '"') { value += '"'; i++ }
      else quoted = !quoted
    } else if (!quoted && c === '#') break
    else if (!quoted && c === ',') { fields.push(value); value = '' }
    else value += c
  }
  expect(quoted).toBe(false)
  fields.push(value)
  return fields
}

function labels(count: number, multiline = true): string[] {
  let state = 23
  const alphabet = ['a', ',', '"', '#', ' ', 'é', '菌', '🧬', '\t', ...(multiline ? ['\r', '\n'] : [])]
  const next = () => (state = (Math.imul(state, 1664525) + 1013904223) >>> 0)
  return Array.from({ length: count }, (_, i) => {
    const length = next() % 40
    return `taxon_${i}:` + Array.from({ length }, () => alphabet[next() % alphabet.length]).join('')
  })
}

test('seeded CSV records preserve delimiters, comments, Unicode and line breaks', () => {
  const values = labels(1024)
  for (let i = 0; i < values.length; i += 4) {
    const fields = values.slice(i, i + 4)
    expect(readRecord(fields.map(csvField).join(','))).toEqual(fields)
  }
  expect(readRecord([csvField('#first'), csvField('middle#tag'), csvField('last')].join(',')))
    .toEqual(['#first', 'middle#tag', 'last'])
})

const rscript = Bun.which('Rscript')
test.skipIf(!rscript)('documented R reader preserves CSV identifiers with comment markers', () => {
  const dir = mkdtempSync(join(tmpdir(), 'metamanifold-csv-'))
  const csv = join(dir, 'results.csv'), expected = join(dir, 'expected.txt')
  const values = ['#first', 'middle#tag', ...labels(128, false)]
  try {
    writeFileSync(csv, '# Differential abundance\r\ntaxon\r\n' + values.map(csvField).join('\r\n') + '\r\n')
    writeFileSync(expected, values.join('\n') + '\n')
    const child = Bun.spawnSync([rscript!, '-e',
      'a<-commandArgs(TRUE); x<-read.csv(a[1],comment.char="#",stringsAsFactors=FALSE); expected<-readLines(a[2],encoding="UTF-8"); stopifnot(identical(enc2utf8(x$taxon),expected))',
      csv, expected])
    expect(new TextDecoder().decode(child.stderr)).toBe('')
    expect(child.exitCode).toBe(0)
  } finally {
    unlinkSync(csv); unlinkSync(expected); rmdirSync(dir)
  }
})
