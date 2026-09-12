// © 2026 Joshua Benjamin Jewell. All rights reserved.
// Licensed under the GNU Affero General Public License version 3 (AGPLv3).


export const PLACE_SELECTED = '#c2255c'
export const SELECT_COLOUR = '#228be6'

export const fmt = (x: number, digits = 4) => Number.isFinite(x) ? Number(x.toPrecision(digits)).toString() : ''
