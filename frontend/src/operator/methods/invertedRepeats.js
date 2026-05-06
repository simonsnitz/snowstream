// Port of findImperfectPalindromes + findBestPalindrome in fetch_operator.py.
import { complement } from '../consensus.js'

const DEFAULT_SPACER_PENALTY = {
  0: 4, 1: 4, 2: 4, 3: 4, 4: 4, 5: 2, 6: 2, 7: 0, 8: 0, 9: -2, 10: -2,
  11: -4, 12: -4, 13: -6, 14: -6, 15: -8, 16: -8, 17: -10, 18: -10,
  19: -12, 20: -12,
}

function findImperfectPalindromes(intergenic, size, winScore, lossScore, spacerPenalty) {
  const len = intergenic.length
  const IRs = len + 1 - 2 * size
  const allScores = []
  for (let i = 0; i < IRs; i++) {
    const repeat = intergenic.slice(i, i + size)
    for (let j = 0; j < IRs; j++) {
      if (j >= 21) continue
      const otherStart = i + size + j
      const otherEnd = i + 2 * size + j
      if (otherEnd > len) continue
      const compare = complement(intergenic.slice(otherStart, otherEnd))
        .split('')
        .reverse()
        .join('')
      if (compare.length !== size) continue
      let score = 0
      for (let k = 0; k < size; k++) {
        score += repeat[k] === compare[k] ? winScore : lossScore
      }
      score += spacerPenalty[String(j)] ?? 0
      const seq =
        repeat +
        intergenic.slice(i + size, i + size + j).toLowerCase() +
        complement(compare).split('').reverse().join('')
      allScores.push({ seq, score })
    }
  }
  if (allScores.length === 0) return null
  const best = Math.max(...allScores.map((s) => s.score))
  return allScores.filter((s) => s.score === best)
}

function findBestPalindromes(intergenic, shortest, longest, winScore, lossScore, spacerPenalty) {
  const upper = intergenic.toUpperCase()
  const ops = []
  for (let i = shortest; i < longest; i++) {
    const r = findImperfectPalindromes(upper, i, winScore, lossScore, spacerPenalty)
    if (r) ops.push(...r)
  }
  if (ops.length === 0) return []
  const max = Math.max(...ops.map((o) => o.score))
  return ops.filter((o) => o.score === max)
}

export const params = {
  win_score: 2,
  loss_score: -2,
  min_operator_length: 5,
  max_operator_length: 15,
  spacer_penalty: DEFAULT_SPACER_PENALTY,
}

// Returns candidate operator sequences for the operator-extraction pipeline.
export function run(promoters, p) {
  const reference = promoters.find((s) => s != null)
  if (!reference) return []
  return findBestPalindromes(
    reference,
    p.min_operator_length,
    p.max_operator_length,
    p.win_score,
    p.loss_score,
    p.spacer_penalty,
  )
}

export const label = 'Inverted repeats'
