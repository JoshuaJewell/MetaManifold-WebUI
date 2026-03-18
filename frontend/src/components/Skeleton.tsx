/** A simple skeleton loading placeholder. */
export function Skeleton({ lines = 3 }: { lines?: number }) {
  const widths = ['60%', '80%', '45%', '70%', '55%']
  return (
    <div className="skeleton">
      {Array.from({ length: lines }, (_, i) => (
        <div key={i} className="skeleton-line" style={{ width: widths[i % widths.length] }} />
      ))}
    </div>
  )
}
