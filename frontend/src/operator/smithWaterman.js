// Local pairwise alignment with affine gap penalties (Gotoh's variant of
// Smith-Waterman). Mirrors Bio.pairwise2.align.localms convention: match
// and mismatch are signed scores; gapOpen and gapExtend are signed
// penalties (typically negative). Opening a gap of length k costs
// gapOpen + k * gapExtend.
//
// Returns the single best local alignment as `{ alignedA, alignedB, score, start, end }`,
// where start/end are the 0-indexed slice bounds in the *original* `a`
// such that a.slice(start, end) covers the aligned region (excluding leading
// gap-only columns in `a`).

const NEG_INF = -1e18

export function localAlign(a, b, { match, mismatch, gapOpen, gapExtend }) {
  const n = a.length
  const m = b.length
  if (n === 0 || m === 0) return null

  // M[i][j]: best score ending at (i,j) with a[i-1] aligned to b[j-1].
  // X[i][j]: best score ending at (i,j) with a[i-1] aligned to a gap (gap in b).
  // Y[i][j]: best score ending at (i,j) with b[j-1] aligned to a gap (gap in a).
  const M = Array.from({ length: n + 1 }, () => new Float64Array(m + 1))
  const X = Array.from({ length: n + 1 }, () => new Float64Array(m + 1))
  const Y = Array.from({ length: n + 1 }, () => new Float64Array(m + 1))
  // Traceback pointer per matrix: 0=stop/zero, 1=M, 2=X, 3=Y.
  const TM = Array.from({ length: n + 1 }, () => new Uint8Array(m + 1))
  const TX = Array.from({ length: n + 1 }, () => new Uint8Array(m + 1))
  const TY = Array.from({ length: n + 1 }, () => new Uint8Array(m + 1))

  for (let i = 0; i <= n; i++) {
    X[i][0] = NEG_INF
    Y[i][0] = NEG_INF
  }
  for (let j = 0; j <= m; j++) {
    X[0][j] = NEG_INF
    Y[0][j] = NEG_INF
  }

  let bestScore = 0
  let bestI = 0
  let bestJ = 0
  let bestMatrix = 1 // M

  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      const s = a[i - 1] === b[j - 1] ? match : mismatch

      // X: gap in b (extend a). Bio.pairwise2 convention: opening a gap of
      // length 1 costs gapOpen only; each additional position costs gapExtend.
      const xOpen = M[i - 1][j] + gapOpen
      const xExt = X[i - 1][j] + gapExtend
      if (xOpen >= xExt) {
        X[i][j] = xOpen
        TX[i][j] = 1 // came from M
      } else {
        X[i][j] = xExt
        TX[i][j] = 2 // came from X
      }

      // Y: gap in a (extend b).
      const yOpen = M[i][j - 1] + gapOpen
      const yExt = Y[i][j - 1] + gapExtend
      if (yOpen >= yExt) {
        Y[i][j] = yOpen
        TY[i][j] = 1 // came from M
      } else {
        Y[i][j] = yExt
        TY[i][j] = 3 // came from Y
      }

      // M: match/mismatch from any predecessor matrix, or restart (0).
      const mFromM = M[i - 1][j - 1] + s
      const mFromX = X[i - 1][j - 1] + s
      const mFromY = Y[i - 1][j - 1] + s
      let bestPred = 0
      let predOrigin = 0
      if (mFromM > bestPred) {
        bestPred = mFromM
        predOrigin = 1
      }
      if (mFromX > bestPred) {
        bestPred = mFromX
        predOrigin = 2
      }
      if (mFromY > bestPred) {
        bestPred = mFromY
        predOrigin = 3
      }
      M[i][j] = bestPred
      TM[i][j] = predOrigin

      if (M[i][j] > bestScore) {
        bestScore = M[i][j]
        bestI = i
        bestJ = j
        bestMatrix = 1
      }
    }
  }

  if (bestScore <= 0) return null

  // Traceback from (bestI, bestJ) in M.
  let i = bestI
  let j = bestJ
  let mat = bestMatrix
  const aOut = []
  const bOut = []
  const endA = bestI

  while (i > 0 && j > 0) {
    if (mat === 1) {
      // M cell: a[i-1] aligned to b[j-1]. Always emit this pair, then
      // either follow the predecessor matrix or stop if origin === 0
      // (i.e., this cell IS the start of the local alignment).
      const origin = TM[i][j]
      aOut.push(a[i - 1])
      bOut.push(b[j - 1])
      i--
      j--
      if (origin === 0) break
      mat = origin
    } else if (mat === 2) {
      // X cell: gap in b
      const origin = TX[i][j]
      aOut.push(a[i - 1])
      bOut.push('-')
      i--
      mat = origin
    } else if (mat === 3) {
      // Y cell: gap in a
      const origin = TY[i][j]
      aOut.push('-')
      bOut.push(b[j - 1])
      j--
      mat = origin
    } else {
      break
    }
  }

  // First indices in raw `a` and `b` that the local alignment matched.
  const startA = i
  const startB = j

  return {
    alignedA: aOut.reverse().join(''),
    alignedB: bOut.reverse().join(''),
    score: bestScore,
    start: startA,
    end: endA,
    startB,
  }
}
