const UNITS: [number, string, string][] = [
  [60,       'second', 'seconds'],
  [60,       'minute', 'minutes'],
  [24,       'hour',   'hours'],
  [30,       'day',    'days'],
  [12,       'month',  'months'],
  [Infinity, 'year',   'years'],
]

/** Parses a server timestamp. The server sends UTC, some without a zone suffix. */
export function parseUtc(dateStr: string): Date {
  return new Date(/[zZ]$|[+-]\d\d:?\d\d$/.test(dateStr) ? dateStr : dateStr + 'Z')
}

/** Returns a human-readable relative time string like "2 min ago". */
export function timeAgo(dateStr: string | null | undefined): string {
  if (!dateStr) return ''
  const date = parseUtc(dateStr)
  if (isNaN(date.getTime())) return dateStr

  let seconds = Math.floor((Date.now() - date.getTime()) / 1000)
  if (seconds < 5) return 'just now'

  for (const [divisor, singular, plural] of UNITS) {
    if (seconds < divisor) {
      const n = Math.floor(seconds)
      return `${n} ${n === 1 ? singular : plural} ago`
    }
    seconds /= divisor
  }
  return dateStr
}
