// SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
// SPDX-License-Identifier: MPL-2.0
import { describe, expect, test } from 'bun:test'
import { exactIntegerText, proportionText } from './validatedStatistics'

describe('descriptive numbers survive the browser boundary', () => {
  test('counts above 2^53 remain exact text', () => {
    expect(exactIntegerText('9007199254740993')).toBe('9007199254740993')
    expect(exactIntegerText('170141183460469231731687303715884105729')).toBe('170141183460469231731687303715884105729')
    expect(exactIntegerText(9007199254740991)).toBe('9007199254740991')
    expect(() => exactIntegerText(9007199254740992)).toThrow()
    for (const bad of ['1e30', '1.0', '-1', 'NaN', '']) expect(() => exactIntegerText(bad)).toThrow()
  })
  test('thirds and zero totals are never rounded into misleading values', () => {
    expect(proportionText('1/3')).toBe('1/3')
    expect(proportionText('0/1')).toBe('0/1')
    expect(proportionText(null)).toBe('Undefined (zero total)')
    expect(() => proportionText('1/0')).toThrow()
    expect(() => proportionText('0.333333')).toThrow()
  })
})
