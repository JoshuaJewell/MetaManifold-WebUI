export function CardActions({ onRename, onDelete }: { onRename: () => void; onDelete: () => void }) {
  return (
    <div style={{
      display: 'flex', gap: 6, justifyContent: 'flex-end',
      marginTop: 10, paddingTop: 8,
      borderTop: '1px solid var(--color-border-light)',
    }}>
      <button
        className="btn"
        style={{ padding: '2px 8px', fontSize: '.78rem' }}
        onClick={e => { e.preventDefault(); e.stopPropagation(); onRename() }}
      >
        Rename
      </button>
      <button
        className="btn"
        style={{ padding: '2px 8px', fontSize: '.78rem', color: '#c92a2a', borderColor: '#ffc9c9' }}
        onClick={e => { e.preventDefault(); e.stopPropagation(); onDelete() }}
      >
        Delete
      </button>
    </div>
  )
}
