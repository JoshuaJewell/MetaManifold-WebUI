const UNITS: [number, string, string][] = [
  [60,       'second', 'seconds'],
  [60,       'minute', 'minutes'],
  [24,       'hour',   'hours'],
  [30,       'day',    'days'],
  [12,       'month',  'months'],
  [Infinity, 'year',   'years'],
]

/** Returns a human-readable relative time string like "2 min ago". */
export function timeAgo(dateStr: string | null | undefined): string {
  if (!dateStr) return ''
  const date = new Date(dateStr)
  if (isNaN(date.getTime())) return dateStr

  let seconds = Math.floor((Date.now() - date.getTime()) / 1000)
  if (seconds < 5) return 'just now'
  if (seconds < 0) return 'just now'

  for (const [divisor, singular, plural] of UNITS) {
    if (seconds < divisor) {
      const n = Math.floor(seconds)
      return `${n} ${n === 1 ? singular : plural} ago`
    }
    seconds /= divisor
  }
  return dateStr
}
