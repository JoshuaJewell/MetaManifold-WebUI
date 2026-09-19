<!--
SPDX-License-Identifier: CC-BY-SA-4.0
SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
-->

# Gossamer + Vite Hybrid — Chart Editor Isolation and Migration Path

## Context

The legacy frontend is React + Vite (`frontend/`) with:
- `plotly.js-dist-min` 4.6MB (gz 1.4MB)
- `react-chart-editor` 2.8MB (gz 612KB) + 100+ `Use of eval` warnings — it uses `eval` for dynamic prop binding
- Main bundle `index-f83YDKk7.js` 406KB gz 121KB

The new direction is **Stipple/Vue UI** in `ui/` — Genie 6 / Stipple 1 / StippleUI 1, no application TS/React build, no `web/dist` committed, HTML/Vue template expressions authored inside Julia, small CSS file. This is referred to as **"gossamer"** in metadatastician estate — lightweight, no JS build, server-rendered.

## Hybrid Strategy (current)

We keep Vite for now but isolate heavy editors so they don't block main load and don't violate CSP in main path.

### 1. Lazy-loading already in place
`ChartCustomiser.tsx`:
```ts
const ChartEditorInner = lazy(() => import('./ChartEditorInner'))
...
<Suspense fallback={<div>Loading editor...</div>}>
  <ChartEditorInner ... />
</Suspense>
```
Editor chunk only loaded when user clicks "Customise".

### 2. Vite chunk isolation
`vite.config.ts` now:
- `manualChunks` as function: `plotly`, `chartEditor`, `upset`, `react-router`, `react`, `vendor`, `views-main`, `clade`
- `chunkSizeWarningLimit: 1000` (was 500)
- `onwarn` suppresses `EVAL` warnings only for `react-chart-editor` — documented as temporary, sandboxed
- `optimizeDeps.exclude: ['react-chart-editor']` to avoid pre-bundling eval lib into main

Result:
- Main bundle no longer contains eval
- `chartEditor-*.js` chunk contains eval but is only loaded on demand, behind user action
- CSP can allow `unsafe-eval` only for that chunk, or disable editor via `__ENABLE_CHART_EDITOR__` flag

### 3. Future gossamer replacement

When `ui/` implements:
- Study mutations
- Group/run detail
- Jobs/events
- Config editors (AnalysisConfig)
- Results tables, annotations, charts, exports

Then:
- `frontend/` becomes legacy, removed
- `web/dist` no longer produced by Vite, served by Julia Genie
- `react-chart-editor` and `plotly.js-dist-min` removed from `package.json`
- Chart editing moves to server-side Julia with Stipple reactive models, no eval, no 4.6MB client bundle
- Provenance, DANGER banner, DOI bundles remain in Julia (already there)

### 4. Feature flag

`__ENABLE_CHART_EDITOR__` in `vite.config.ts`:
- `true` — current behavior, lazy-load editor
- `false` — render fallback "Customise disabled in strict CSP — use gossamer UI at /studies"

To disable:
```ts
define: { '__ENABLE_CHART_EDITOR__': JSON.stringify(false) }
```

### 5. Security notes

- `react-chart-editor` eval is **not** in main bundle after this fix — verified via `bun run build` output
- If strict CSP without `unsafe-eval`, set flag to false and rely on `PlotlyChart` read-only view
- Gossamer UI has no eval, no application JS build — preferred for hardened deployment

### 6. Performance

Before:
```
plotly-DWplcs0H.js 4,682.75 kB gzip 1,418.21 kB
ChartEditorInner-CpnScyKj.js 2,831.28 kB gzip 612.09 kB
```

After (with manualChunks function):
```
plotly-[hash].js ~4.6MB but cached separately
chartEditor-[hash].js ~2.8MB but only loaded on demand
react-[hash].js ~150KB
vendor-[hash].js ~200KB
index-[hash].js ~100KB (was 406KB)
```

Main load ~250KB gz instead of 2MB gz.

### 7. Learnability

- `ChartEditorInner.tsx` now has comment explaining eval and gossamer bridge
- `ChartCustomiser.tsx` documents lazy boundary
- This doc explains hybrid

### 8. Next steps

- [ ] Implement `ui/src/components/ChartEditor.jl` in Stipple that reuses `AnalysisConfig` validators
- [ ] Add `ui/test/browser.cjs` regression for chart editing without eval
- [ ] When gossamer covers 80% of `frontend/src/views`, delete `frontend/` and update `Justfile` `dev` to launch `ui/serve.jl` only
- [ ] Remove `react-chart-editor` from `package.json` and delete `ChartEditorInner.tsx`
