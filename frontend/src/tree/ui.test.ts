// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

import { describe, expect, test } from 'bun:test'
import { taxonRuns } from './ui'

const italics = (label: string) => taxonRuns(label).map(r => (r.italic ? `*${r.text}*` : r.text)).join('')

describe('taxonRuns', () => {
  test('genus and species are italic, accessions and subtypes upright', () => {
    expect(italics('Myocastor coypus PZ941397')).toBe('*Myocastor coypus* PZ941397')
    expect(italics('G. intestinalis AII KM926505')).toBe('*G. intestinalis* AII KM926505')
    expect(italics('ST5 MK375236')).toBe('ST5 MK375236')
  })
  test('qualifiers and joined names', () => {
    expect(italics('Periplaneta americana isolates')).toBe('*Periplaneta americana* isolates')
    expect(italics('Retortamonas sp. AF439348')).toBe('*Retortamonas* sp. AF439348')
    expect(italics('Holomastigotoides + Spirotrichonymphida')).toBe('*Holomastigotoides* + *Spirotrichonymphida*')
    expect(italics('Castor fiber')).toBe('*Castor fiber*')
  })
  test('configured exceptions stay upright', () => {
    expect(taxonRuns('Outgroup', new Set(['Outgroup']))).toEqual([{ text: 'Outgroup', italic: false }])
    expect(italics('Outgroup')).toBe('*Outgroup*')
  })
  test('the runs rebuild the label', () => {
    const label = 'Platysternon megacephalum megacephalum KT438715'
    expect(taxonRuns(label).map(r => r.text).join('')).toBe(label)
  })
})
