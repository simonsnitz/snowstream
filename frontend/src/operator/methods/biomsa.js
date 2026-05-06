// BioMSA-based candidate operator finder.
//
// Strategy: align all homolog promoters with biomsa, then for each window
// length L in [min_operator_length, max_operator_length] scan the alignment
// and return the consensus of the most-conserved L-column window. The
// shared operator pipeline then SW-aligns each consensus to every homolog
// promoter and picks the one yielding the best overall consensus_score.

import biomsa from 'biomsa'

export const params = {
  min_operator_length: 15,
  max_operator_length: 25,
}

export const label = 'BioMSA'

// Per-column consensus + conservation score over an MSA.
// score = fraction of non-gap rows that match the most-common non-gap base.
function columnStats(msa, col) {
  const counts = Object.create(null)
  let nonGap = 0
  for (let r = 0; r < msa.length; r++) {
    const ch = msa[r][col]
    if (ch === '-') continue
    counts[ch] = (counts[ch] ?? 0) + 1
    nonGap += 1
  }
  if (nonGap === 0) return { base: '-', score: 0 }
  let bestBase = '-'
  let bestCount = 0
  for (const [b, c] of Object.entries(counts)) {
    if (c > bestCount) {
      bestCount = c
      bestBase = b
    }
  }
  return { base: bestBase, score: bestCount / nonGap }
}

// Best window of a given length: returns { start, mean_score, consensus }
// or null if the alignment is shorter than the window.
function bestWindow(msa, perColumn, length) {
  if (perColumn.length < length) return null
  let bestStart = 0
  let bestSum = -1
  let runningSum = 0
  for (let i = 0; i < length; i++) runningSum += perColumn[i].score
  bestSum = runningSum
  for (let start = 1; start + length <= perColumn.length; start++) {
    runningSum += perColumn[start + length - 1].score - perColumn[start - 1].score
    if (runningSum > bestSum) {
      bestSum = runningSum
      bestStart = start
    }
  }
  // Consensus over the chosen window, dropping any all-gap columns.
  let consensus = ''
  for (let i = bestStart; i < bestStart + length; i++) {
    if (perColumn[i].base !== '-') consensus += perColumn[i].base
  }
  return { start: bestStart, mean_score: bestSum / length, consensus }
}

export async function run(promoters, p = {}) {
  const seqs = promoters.filter((s) => s != null && s.length > 0)
  if (seqs.length < 2) return []

  const minL = p.min_operator_length ?? params.min_operator_length
  const maxL = p.max_operator_length ?? params.max_operator_length

  let aligned
  try {
    aligned = await biomsa.align(seqs, { type: 'nucleic' })
  } catch {
    return []
  }
  if (!aligned?.length) return []

  const cols = aligned[0].length
  const perColumn = new Array(cols)
  for (let c = 0; c < cols; c++) perColumn[c] = columnStats(aligned, c)

  const seenSeqs = new Set()
  const candidates = []
  for (let L = minL; L <= maxL; L++) {
    const w = bestWindow(aligned, perColumn, L)
    if (!w || !w.consensus) continue
    if (seenSeqs.has(w.consensus)) continue
    seenSeqs.add(w.consensus)
    candidates.push({ seq: w.consensus })
  }
  return candidates
}
