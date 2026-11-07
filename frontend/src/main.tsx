// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
// SPDX-License-Identifier: MPL-2.0

import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import { loadConfig } from './api/client'
import { App } from './App'
import './styles/app.css'

loadConfig().then(() => {
  createRoot(document.getElementById('root')!).render(
    <StrictMode>
      <App />
    </StrictMode>,
  )
})
