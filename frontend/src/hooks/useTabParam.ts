// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

import { useCallback } from 'react'
import { useSearchParams } from 'react-router-dom'

/**
 * A tab selection kept in the URL query (`?key=value`), so it survives reloads
 * and can be linked to. Unknown values fall back to `fallback`.
 */
export function useTabParam<T extends string>(key: string, values: readonly T[], fallback: T) {
  const [params, setParams] = useSearchParams()
  const raw = params.get(key)
  const value = (values as readonly string[]).includes(raw ?? '') ? (raw as T) : fallback
  const setValue = useCallback((next: T) => {
    setParams(prev => {
      const p = new URLSearchParams(prev)
      if (next === fallback) p.delete(key)
      else p.set(key, next)
      return p
    }, { replace: true })
  }, [key, fallback, setParams])
  return [value, setValue] as const
}
