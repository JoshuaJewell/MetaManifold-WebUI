// SPDX-License-Identifier: AGPL-3.0-only
// SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell; 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
// Public surface of the exploratory dendrogram-receipts extension.
export * from './types'
export { parseBridgeIndex, parseReceipt, resolveSelection, BridgeContractError, type Selection } from './bridgeIndex'
export { attachDendrogram, createReceiptController, type PanelState } from './controller'
export { renderPanel } from './renderPanel'
