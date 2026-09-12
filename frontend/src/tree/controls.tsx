// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).
import type { NodeStyle, Rgba } from './model'
import type { SymbolStyle } from './layout'
import styles from './Tree.module.css'

export function RgbaInput({ value, onChange, label }: { value: Rgba; onChange: (v: Rgba) => void; label: string }) {
  return (
    <span className={styles.rgba} title={`${label}: colour and opacity`}>
      <input type="color" aria-label={`${label} colour`} value={value.hex} onChange={e => onChange({ ...value, hex: e.target.value })} />
      <input type="number" aria-label={`${label} opacity`} min={0} max={1} step={0.05} value={value.alpha} className={styles.alpha}
        onChange={e => onChange({ ...value, alpha: Math.min(1, Math.max(0, Number(e.target.value))) })} />
    </span>
  )
}

export function DatasetControls({ title, style, onChange, children }: {
  title: string; style: SymbolStyle; onChange: (p: Partial<SymbolStyle>) => void; children?: React.ReactNode
}) {
  return (
    <div className={styles.dataset}>
      <div className={styles.datasetHead}>
        <strong>{title}</strong>
        <label><input type="checkbox" checked={style.show} onChange={e => onChange({ show: e.target.checked })} /> Symbols</label>
      </div>
      {style.show && (
        <div className={styles.datasetGrid}>
          <span>Shape</span>
          <select value={style.shape} onChange={e => onChange({ shape: e.target.value as SymbolStyle['shape'] })}>
            <option value="circle">Circle</option>
            <option value="square">Square</option>
            <option value="diamond">Diamond</option>
            <option value="triangle">Triangle</option>
          </select>
          <span>Fill</span>
          <RgbaInput label="Fill" value={style.fill} onChange={fill => onChange({ fill })} />
          <span>Outline</span>
          <RgbaInput label="Outline" value={style.stroke} onChange={stroke => onChange({ stroke })} />
          <span>Size (px)</span>
          <span className={styles.group}>
            <input type="number" aria-label="Smallest size" min={0} max={100} value={style.minSize} className={styles.num}
              onChange={e => onChange({ minSize: Math.max(0, Number(e.target.value) || 0) })} />
            to
            <input type="number" aria-label="Largest size" min={1} max={200} value={style.maxSize} className={styles.num}
              onChange={e => onChange({ maxSize: Math.max(1, Number(e.target.value) || 1) })} />
          </span>
        </div>
      )}
      <div className={styles.datasetExtra}>{children}</div>
    </div>
  )
}

export function StyleEditor({ style, onChange }: { style: NodeStyle | null; onChange: (s: NodeStyle | null) => void }) {
  const cur: NodeStyle = style ?? { inherit: true }
  const update = (patch: Partial<NodeStyle>) => {
    const next: NodeStyle = { ...cur, ...patch }
    for (const k of ['bold', 'italic', 'label', 'branch', 'fill'] as const) if (next[k] === undefined) delete next[k]
    const empty = next.bold === undefined && next.italic === undefined && !next.label && !next.branch && !next.fill
    onChange(empty ? null : next)
  }
  const tri = (v: boolean | undefined) => v === undefined ? '' : v ? 'on' : 'off'
  const fromTri = (v: string) => v === '' ? undefined : v === 'on'
  return (
    <div className={styles.datasetGrid}>
      <span>Bold</span>
      <select value={tri(cur.bold)} onChange={e => update({ bold: fromTri(e.target.value) })}>
        <option value="">inherit</option><option value="on">on</option><option value="off">off</option>
      </select>
      <span>Italic</span>
      <select value={tri(cur.italic)} onChange={e => update({ italic: fromTri(e.target.value) })}>
        <option value="">inherit</option><option value="on">on</option><option value="off">off</option>
      </select>
      <span>Label colour</span>
      <span className={styles.group}>
        <input type="checkbox" aria-label="Set label colour" checked={!!cur.label}
          onChange={e => update({ label: e.target.checked ? { hex: '#1971c2', alpha: 1 } : undefined })} />
        {cur.label && <RgbaInput label="Label colour" value={cur.label} onChange={label => update({ label })} />}
      </span>
      <span>Branch colour</span>
      <span className={styles.group}>
        <input type="checkbox" aria-label="Set branch colour" checked={!!cur.branch}
          onChange={e => update({ branch: e.target.checked ? { hex: '#c92a2a', alpha: 1 } : undefined })} />
        {cur.branch && <RgbaInput label="Branch colour" value={cur.branch} onChange={branch => update({ branch })} />}
      </span>
      <span>Triangle fill</span>
      <span className={styles.group}>
        <input type="checkbox" aria-label="Set triangle fill" checked={!!cur.fill}
          onChange={e => update({ fill: e.target.checked ? { hex: cur.branch?.hex ?? '#868e96', alpha: 0.3 } : undefined })} />
        {cur.fill && <RgbaInput label="Triangle fill" value={cur.fill} onChange={fill => update({ fill })} />}
      </span>
      <span>Descendants</span>
      <label><input type="checkbox" checked={cur.inherit} onChange={e => update({ inherit: e.target.checked })} /> inherit this style</label>
    </div>
  )
}
