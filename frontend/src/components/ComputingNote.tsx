// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import { useEffect, useState } from 'react'

// Shown under an analysis's buttons while its request runs. The R-backed
// analyses can take tens of seconds on a large study, and a disabled button
// alone reads as a page that has stopped responding.
export function ComputingNote({ active, label = 'Computing' }: { active: boolean; label?: string }) {
  const [seconds, setSeconds] = useState(0)

  useEffect(() => {
    if (!active) { setSeconds(0); return }
    const started = Date.now()
    const timer = setInterval(() => setSeconds(Math.floor((Date.now() - started) / 1000)), 1000)
    return () => clearInterval(timer)
  }, [active])

  if (!active) return null
  return (
    <p className="loading computing-note" role="status">
      <span className="computing-spinner" aria-hidden="true" />
      {label}…{seconds >= 2 && ` ${seconds}s`}
      {seconds >= 10 && ' (large studies can take a minute or more)'}
    </p>
  )
}
