// SPDX-License-Identifier: AGPL-3.0-only
import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// Hybrid Vite + Gossamer strategy:
// - Current React/Vite app remains default, but heavy editors are isolated via lazy() + manualChunks
// - react-chart-editor uses eval internally (100+ warnings) — it is sandboxed in its own chunk
//   and only loaded when user clicks "Customise". This is a temporary bridge until
//   the Stipple/Vue UI in ui/ (Genie/Stipple/StippleUI, no JS build) — referred to as
//   "gossamer" in metadatastician estate — replaces the React chart editing path.
// - Plotly is also isolated (4.6MB) and lazy-loadable in future; for now manualChunk keeps it out of main.
// - When gossamer is ready, this Vite setup becomes legacy and can be removed; web/dist will be served by Julia.
// See docs/migration/gossamer-vite-hybrid.md for migration plan.

export default defineConfig({
  plugins: [react()],
  define: {
    // Feature flag: disable heavy editor in strict CSP environments
    '__ENABLE_CHART_EDITOR__': JSON.stringify(true),
  },
  build: {
    outDir:      '../web/dist',
    emptyOutDir: true,
    sourcemap:   false,
    chunkSizeWarningLimit: 1000,
    rollupOptions: {
      // Suppress eval warnings only for react-chart-editor — known, isolated, temporary
      onwarn(warning, warn) {
        // react-chart-editor lib uses eval for dynamic prop binding — sandboxed in its own chunk
        // See https://github.com/plotly/react-chart-editor/issues
        if (warning.code === 'EVAL' && warning.id && warning.id.includes('react-chart-editor')) {
          return
        }
        // Also suppress circular dependency warnings from plotly
        if (warning.code === 'CIRCULAR_DEPENDENCY' && warning.ids?.some((id: string) => id.includes('plotly'))) {
          return
        }
        warn(warning)
      },
      output: {
        // Manual chunks for better caching and to isolate heavy deps
        // Each chunk can be independently cached and only loaded when needed
        manualChunks(id) {
          if (id.includes('node_modules')) {
            if (id.includes('plotly.js-dist-min') || id.includes('plotly.js')) {
              return 'plotly'
            }
            if (id.includes('react-chart-editor')) {
              return 'chartEditor'
            }
            if (id.includes('@upsetjs')) {
              return 'upset'
            }
            if (id.includes('react-router')) {
              return 'react-router'
            }
            if (id.includes('react-dom') || id.includes('react/')) {
              return 'react'
            }
            // Other vendor libs
            return 'vendor'
          }
          // App code splitting: isolate views that are heavy
          if (id.includes('src/views/RunView') || id.includes('src/views/StudyView')) {
            return 'views-main'
          }
          if (id.includes('src/components/CladeCumulus')) {
            return 'clade'
          }
        },
        // Ensure eval-containing chunk is marked as having dynamic requires
        // and is not inlined into main
        chunkFileNames: 'assets/[name]-[hash].js',
        entryFileNames: 'assets/[name]-[hash].js',
      },
    },
  },
  optimizeDeps: {
    // Exclude heavy deps that are lazy-loaded to avoid pre-bundling them into main
    exclude: ['react-chart-editor'],
    include: ['react', 'react-dom', 'react-router-dom', 'plotly.js-dist-min'],
  },
  server: {
    proxy: {
      // SSE endpoint needs special handling to prevent buffering
      '/api/v1/events': {
        target: 'http://127.0.0.1:8080',
        // Disable response buffering so SSE frames stream through immediately
        configure: (proxy) => {
          proxy.on('proxyRes', (proxyRes) => {
            proxyRes.headers['cache-control'] = 'no-cache'
            proxyRes.headers['x-accel-buffering'] = 'no'
          })
        },
      },
      '/api':   'http://127.0.0.1:8080',
      '/files': 'http://127.0.0.1:8080',
    },
  },
})
