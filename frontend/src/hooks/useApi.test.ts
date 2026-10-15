import { describe, expect, test } from 'bun:test'
import { failedState, type State } from './useApi'

const shown: State<string[]> = { data: ['row from the old filters'], loading: true, error: null }

describe('useApi failed request', () => {
  test('drops data fetched for other inputs', () => {
    expect(failedState(shown, 'HTTP 500', false)).toEqual({ data: null, loading: false, error: 'HTTP 500' })
  })

  test('keeps data when the same request is re-run and fails', () => {
    expect(failedState(shown, 'HTTP 500', true)).toEqual({
      data: ['row from the old filters'], loading: false, error: 'HTTP 500',
    })
  })
})
