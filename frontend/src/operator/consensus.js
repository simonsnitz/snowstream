// Ports of getConsensus, get_consensus_score, generate_frequency_matrix
// from src/fetch_operator.py.

import { localAlign } from './smithWaterman.js'

export function complement(seq) {
  const map = { A: 'T', T: 'A', C: 'G', G: 'C' }
  let out = ''
  for (const c of seq.toUpperCase()) {
    if (!(c in map)) break
    out += map[c]
  }
  return out
}

export function reverseComplement(seq) {
  return complement(seq).split('').reverse().join('')
}

// Locate `operator` within `intergenic` via local alignment, returning the
// extracted operator (with `extLength` lowercase flanking bases on each side).
// Mirrors findOperatorInIntergenic in fetch_operator.py.
export function findOperatorInIntergenic(intergenic, operator, params) {
  const { extension_length, align_match, align_mismatch, gap_open, gap_extend } = params
  const opLen = operator.length

  const result = localAlign(intergenic, operator, {
    match: align_match,
    mismatch: align_mismatch,
    gapOpen: gap_open,
    gapExtend: gap_extend,
  })
  if (!result) return null

  const maxScore = 2 * opLen
  const cutoff = maxScore * 0.1
  if (result.score <= cutoff) return null

  // Anchor the operator at its full length: even if the local match skipped
  // some leading chars of `operator` (because they didn't score positively),
  // Bio.pairwise2 reports the operator content extending into those leading
  // columns of intergenic. We mirror that behaviour by shifting the start
  // back by `startB` (operator chars before the matched region).
  const absoluteBegin = result.start - result.startB
  const absoluteEnd = absoluteBegin + opLen

  let extracted
  try {
    const upstream = intergenic.slice(absoluteBegin - extension_length, absoluteBegin).toLowerCase()
    const mid = intergenic.slice(absoluteBegin, absoluteEnd)
    const downstream = intergenic.slice(absoluteEnd, absoluteEnd + extension_length).toLowerCase()
    extracted = upstream + mid + downstream
  } catch {
    extracted = intergenic.slice(absoluteBegin, absoluteEnd)
  }

  return { operator: extracted, score: result.score }
}

// Build per-position consensus from a list of equal-length operator strings.
// Returns { motif_data: [{base, score}], num_seqs } where score is the
// fractional conservation at each position (max conservation = 1.0).
export function getConsensus(metrics) {
  const operators = metrics
    .filter((m) => m['Align score'] !== 0)
    .map((m) => m['Predicted operator'])
  if (operators.length === 0) return { motif_data: [], num_seqs: 0 }

  const refLen = operators[0].length
  const baep = Array.from({ length: refLen }, () => Object.create(null))

  for (const op of operators) {
    if (op.length !== refLen) continue
    for (let pos = 0; pos < op.length; pos++) {
      const base = op[pos]
      baep[pos][base] = (baep[pos][base] ?? 0) + 1
    }
  }

  const maxValues = baep.map((counts) => Math.max(...Object.values(counts)))
  const overallMax = Math.max(...maxValues)
  const maxValuesPercent = maxValues.map((v) => Math.round((v / overallMax) * 100) / 100)

  const motif_data = baep.map((counts, i) => {
    const max = maxValues[i]
    const base = Object.entries(counts).find(([, v]) => v === max)?.[0] ?? 'N'
    return { base, score: maxValuesPercent[i] }
  })

  return { motif_data, num_seqs: operators.length }
}

// Quality score for a candidate operator: fraction of "core" (uppercase)
// positions whose conservation^2 sums up to the maximum possible score.
// Mirrors get_consensus_score in fetch_operator.py.
export function getConsensusScore(operator, consensusData, extLength) {
  let maxScore = 0
  let consensusScore = 0
  const motif = consensusData.motif_data
  for (let i = 0; i < operator.length; i++) {
    const ch = operator[i]
    if (ch >= 'A' && ch <= 'Z') {
      maxScore += 1
      const m = motif[i + extLength]
      if (m) consensusScore += m.score ** 2
    }
  }
  if (maxScore === 0) return 0
  return Math.round((consensusScore / maxScore) * 100 * 1000) / 1000
}

// 4xN frequency matrix [A, C, G, T] of equal-length operators in metadata.
// Mirrors generate_frequency_matrix in fetch_operator.py.
export function generateFrequencyMatrix(metadata) {
  const operators = []
  let opLength = null
  for (const h of metadata) {
    const op = (h['Predicted operator'] || '').toUpperCase()
    if (opLength === null && op) opLength = op.length
    if (op.length === opLength && /^[ATCG]+$/.test(op)) operators.push(op)
  }
  if (operators.length === 0 || opLength === null) return []

  const numOps = operators.length
  const matrix = []
  for (let i = 0; i < opLength; i++) {
    const base = [0, 0, 0, 0]
    for (const op of operators) {
      const c = op[i]
      if (c === 'A') base[0] += 1 / numOps
      else if (c === 'C') base[1] += 1 / numOps
      else if (c === 'G') base[2] += 1 / numOps
      else if (c === 'T') base[3] += 1 / numOps
    }
    matrix.push(base.map((x) => Math.round((x + Number.EPSILON) * 100) / 100))
  }
  return matrix
}
