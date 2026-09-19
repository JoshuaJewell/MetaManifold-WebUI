// SPDX-License-Identifier: AGPL-3.0-only
import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  build: {
    outDir:      '../web/dist',
    emptyOutDir: true,
    sourcemap:   false,
    chunkSizeWarningLimit: 1000,
    rollupOptions: {
      output: {
        manualChunks: {
          plotly: ['plotly.js-dist-min'],
          react: ['react', 'react-dom', 'react-router-dom'],
          chartEditor: ['react-chart-editor'],
          upset: ['@upsetjs/react'],
        },
      },
    },
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
