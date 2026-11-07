// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

/** Undo and redo for a document edited in small steps. `at` is when the last edit was
 *  recorded; edits closer together than the merge window, such as the keystrokes of one
 *  typed value, undo as one step. */
export interface History<T> {
  past: T[]
  future: T[]
  at: number
}

export const emptyHistory = <T>(): History<T> => ({ past: [], future: [], at: 0 })

const MERGE_MS = 600
const LIMIT = 100

/** The history after `current` is replaced by an edit made at `now`. */
export function record<T>(h: History<T>, current: T, now: number): History<T> {
  const merge = h.past.length > 0 && now - h.at < MERGE_MS
  const past = merge ? h.past : [...h.past, current].slice(-LIMIT)
  return { past, future: [], at: now }
}

/** The document to show after undoing, and the history that follows; null when there is nothing to undo. */
export function undo<T>(h: History<T>, current: T): { doc: T; history: History<T> } | null {
  if (!h.past.length) return null
  return { doc: h.past[h.past.length - 1], history: { past: h.past.slice(0, -1), future: [current, ...h.future], at: 0 } }
}

/** The document to show after redoing, and the history that follows; null when there is nothing to redo. */
export function redo<T>(h: History<T>, current: T): { doc: T; history: History<T> } | null {
  if (!h.future.length) return null
  return { doc: h.future[0], history: { past: [...h.past, current], future: h.future.slice(1), at: 0 } }
}

/** Ctrl or Cmd with Z undoes; adding Shift, or Ctrl or Cmd with Y, redoes. */
export function undoKey(e: Pick<KeyboardEvent, 'key' | 'ctrlKey' | 'metaKey' | 'shiftKey' | 'altKey'>): 'undo' | 'redo' | null {
  if (!(e.ctrlKey || e.metaKey) || e.altKey) return null
  const k = e.key.toLowerCase()
  if (k === 'z') return e.shiftKey ? 'redo' : 'undo'
  if (k === 'y' && !e.shiftKey) return 'redo'
  return null
}
