// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).

import { useState } from 'react'
import { api } from '../api/client'
import { useToast } from './Toast'
import type { ContaminationStats } from '../api/types'
import { CONTAM_FILTER_FIELDS, FUNCDB_FIELDS } from './annotationShared'

interface ContamFilterState {
  blacklist: Record<string, string>
  whitelist: Record<string, string>
}

export function ContamStatsBar({ stats }: { stats: ContaminationStats }) {
  const fmt = (n: number) => n.toLocaleString()
  const pct = (n: number) =>
    stats.total.rows > 0 ? ((n / stats.total.rows) * 100).toFixed(1) : '0.0'

  const rows: Array<{ key: keyof Omit<ContaminationStats, 'total'>; label: string; color: string }> = [
    { key: 'yes', label: 'Marked contamination', color: '#c92a2a' },
    { key: 'no', label: 'Marked non-contamination', color: '#2b8a3e' },
    { key: 'unassigned', label: 'Unassigned', color: 'var(--color-muted-fg)' },
  ]

  return (
    <div style={{ marginTop: 10, fontSize: '.82rem' }}>
      <div style={{ marginBottom: 6, color: 'var(--color-muted-fg)' }}>
        Total {fmt(stats.total.rows)} rows - {fmt(stats.total.reads)} reads
      </div>
      {rows.map(({ key, label, color }) => {
        const bucket = stats[key]
        const width = parseFloat(pct(bucket.rows))
        return (
          <div key={key} style={{ display: 'flex', alignItems: 'center', gap: 8, marginBottom: 4 }}>
            <div style={{ width: 120, height: 6, background: 'var(--color-border)', borderRadius: 3, flexShrink: 0 }}>
              <div style={{ width: `${width}%`, height: '100%', background: color, borderRadius: 3, transition: 'width .3s' }} />
            </div>
            <span style={{ color, fontWeight: 600, minWidth: 90 }}>{label}</span>
            <span style={{ color: 'var(--color-muted-fg)' }}>
              {fmt(bucket.rows)} rows - {fmt(bucket.reads)} reads - {pct(bucket.rows)}%
            </span>
          </div>
        )
      })}
    </div>
  )
}

export function ContamFilterConfig({
  config,
  onChange,
  onApply,
  applying,
  stats,
}: {
  config: ContamFilterState
  onChange: (config: ContamFilterState) => void
  onApply: () => void
  applying: boolean
  stats: ContaminationStats | null
}) {
  const setField = (side: 'blacklist' | 'whitelist', key: string, value: string) =>
    onChange({ ...config, [side]: { ...config[side], [key]: value } })

  return (
    <div style={{
      marginBottom: 12,
      padding: '10px 12px',
      background: 'var(--color-surface)',
      borderRadius: 6,
      border: '1px solid var(--color-border)',
    }}>
      <div style={{ fontWeight: 600, fontSize: '.82rem', marginBottom: 6 }}>Contamination Filter</div>
      <div style={{ fontSize: '.76rem', color: 'var(--color-muted-fg)', marginBottom: 10 }}>
        Apply filter recomputes contamination from the current contamination and non-contamination filters and overwrites manual row edits.
      </div>
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 16 }}>
        {(['blacklist', 'whitelist'] as const).map(side => (
          <div key={side}>
            <div style={{
              fontWeight: 600,
              fontSize: '.78rem',
              marginBottom: 8,
              color: side === 'blacklist' ? '#c92a2a' : '#2b8a3e',
            }}>
              {side === 'blacklist' ? 'Contamination filter' : 'Non-contamination filter'}
            </div>
            <div style={{ fontSize: '.72rem', color: 'var(--color-muted-fg)', marginBottom: 8 }}>
              {side === 'blacklist' ? 'Matching rows are tagged as contamination.' : 'Matching rows are tagged as non-contamination.'}
            </div>
            {CONTAM_FILTER_FIELDS.map(({ key, label }) => (
              <label key={key} style={{ display: 'flex', flexDirection: 'column', gap: 2, marginBottom: 8, fontSize: '.78rem' }}>
                <span>{label}</span>
                <textarea
                  value={config[side][key] ?? ''}
                  onChange={event => setField(side, key, event.target.value)}
                  rows={3}
                  placeholder="One value per line"
                  style={{
                    resize: 'vertical',
                    padding: '3px 6px',
                    borderRadius: 4,
                    border: '1px solid var(--color-border)',
                    fontSize: '.78rem',
                    background: 'var(--color-bg)',
                    fontFamily: 'monospace',
                  }}
                />
              </label>
            ))}
          </div>
        ))}
      </div>
      <button className="btn btn-primary" onClick={onApply} disabled={applying} style={{ marginTop: 4 }}>
        {applying ? 'Applying...' : 'Apply filter'}
      </button>
      {stats && <ContamStatsBar stats={stats} />}
    </div>
  )
}

export function AddFuncdbModal({
  onClose,
  prefill,
  defaultModifiedBy,
}: {
  onClose: () => void
  prefill?: Record<string, string>
  defaultModifiedBy: string
}) {
  const toast = useToast()
  const [fields, setFields] = useState<Record<string, string>>(prefill ?? {})
  const [modifiedBy, setModifiedBy] = useState(defaultModifiedBy)
  const [saving, setSaving] = useState(false)

  const setField = (key: string, value: string) => setFields(current => ({ ...current, [key]: value }))

  const handleSubmit = async (event: React.FormEvent) => {
    event.preventDefault()
    const hasTaxon = ['Domain', 'supergroup', 'division', 'class', 'order', 'family', 'Genus', 'Species']
      .some(key => fields[key]?.trim())
    if (!hasTaxon) {
      toast.error('At least one taxonomy field is required')
      return
    }
    if (!fields.Function?.trim()) {
      toast.error('Function is required')
      return
    }

    setSaving(true)
    try {
      await api.annotations.addFuncdbEntry(fields, modifiedBy || undefined)
      toast.success('FuncDB entry added')
      onClose()
    } catch (err) {
      toast.error(err instanceof Error ? err.message : 'Failed to add entry')
    } finally {
      setSaving(false)
    }
  }

  return (
    <div
      style={{
        position: 'fixed', inset: 0, zIndex: 1000,
        background: 'rgba(0,0,0,.45)', display: 'flex',
        alignItems: 'center', justifyContent: 'center',
      }}
      onClick={event => { if (event.target === event.currentTarget) onClose() }}
    >
      <div style={{
        background: 'var(--color-bg)', borderRadius: 8, padding: '20px 24px',
        width: 520, maxHeight: '85vh', overflowY: 'auto',
        boxShadow: '0 8px 32px rgba(0,0,0,.25)',
      }}>
        <h3 style={{ margin: '0 0 16px', fontSize: '1rem' }}>Add FuncDB Entry</h3>
        <form onSubmit={handleSubmit}>
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '8px 12px', marginBottom: 12 }}>
            {FUNCDB_FIELDS.map(({ key, label, required }) => (
              <label key={key} style={{ display: 'flex', flexDirection: 'column', gap: 3, fontSize: '.82rem' }}>
                <span style={{ fontWeight: required ? 600 : 400 }}>
                  {label}{required && <span style={{ color: '#c92a2a' }}> *</span>}
                </span>
                <input
                  value={fields[key] ?? ''}
                  onChange={event => setField(key, event.target.value)}
                  style={{
                    padding: '3px 6px',
                    borderRadius: 4,
                    border: '1px solid var(--color-border)',
                    fontSize: '.82rem',
                    background: 'var(--color-bg)',
                  }}
                />
              </label>
            ))}
          </div>
          <label style={{ display: 'flex', flexDirection: 'column', gap: 3, fontSize: '.82rem', marginBottom: 16 }}>
            <span>Modified by</span>
            <input
              value={modifiedBy}
              onChange={event => setModifiedBy(event.target.value)}
              placeholder="Your name"
              style={{
                padding: '3px 6px',
                borderRadius: 4,
                border: '1px solid var(--color-border)',
                fontSize: '.82rem',
                background: 'var(--color-bg)',
                width: '50%',
              }}
            />
          </label>
          <div style={{ display: 'flex', gap: 8, justifyContent: 'flex-end' }}>
            <button type="button" className="btn" onClick={onClose}>Cancel</button>
            <button type="submit" className="btn btn-primary" disabled={saving}>
              {saving ? 'Saving...' : 'Add entry'}
            </button>
          </div>
        </form>
      </div>
    </div>
  )
}
