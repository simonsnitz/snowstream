// Method-agnostic operator-extraction pipeline. Mirrors the body of
// fetch_operator() in src/fetch_operator.py.
//
// Input:
//   homologs: [{ uniprot_id, promoter, ... }]   (promoter may be null)
//   candidates: [{ seq }]                        (from a method's run())
//   params: { extension_length, align_match, align_mismatch, gap_open, gap_extend }
//
// Output:
//   {
//     consensus_score, num_seqs, consensus_seq, native_operator,
//     motif: [{base, score}], aligned_seqs: [{uniprot_id, predicted_operator, align_score}],
//     frequency_matrix: [[a,c,g,t], ...],
//     intergenic
//   }

import {
  findOperatorInIntergenic,
  getConsensus,
  getConsensusScore,
  generateFrequencyMatrix,
} from './consensus.js'

export const DEFAULT_PARAMS = {
  extension_length: 5,
  align_match: 2,
  align_mismatch: -0.5,
  gap_open: -100,
  gap_extend: 0,
}

export function extractOperators(homologs, candidates, params = {}) {
  const p = { ...DEFAULT_PARAMS, ...params }
  const promoters = homologs.map((h) => h.promoter)
  const reference = promoters.find((s) => s != null)
  if (!reference || candidates.length === 0) {
    return null
  }

  const best = {
    consensus_score: 0,
    num_seqs: 0,
    consensus_seq: '',
    native_operator: '',
    motif: [],
    aligned_seqs: [],
    frequency_matrix: [],
    intergenic: reference,
  }

  for (const candidate of candidates) {
    const metrics = []
    for (const h of homologs) {
      if (!h.promoter) continue
      const op = findOperatorInIntergenic(h.promoter, candidate.seq, p)
      if (op) {
        metrics.push({
          'Uniprot Id': h.uniprot_id,
          'Predicted operator': op.operator,
          'Align score': op.score,
        })
      }
    }
    if (metrics.length === 0) continue

    const consensus = getConsensus(metrics)
    const score = getConsensusScore(candidate.seq, consensus, p.extension_length)

    // Pull out the native operator sequence from the original (reference) promoter.
    const nativeOp = findOperatorInIntergenic(reference, candidate.seq, p)
    const nativeSeq = nativeOp ? nativeOp.operator : candidate.seq
    const consensusSeq = consensus.motif_data.map((m) => m.base).join('')

    if (score > best.consensus_score) {
      best.consensus_score = score
      best.native_operator = nativeSeq
      best.consensus_seq = consensusSeq
      best.num_seqs = consensus.num_seqs
      best.motif = consensus.motif_data
      best.aligned_seqs = metrics.map((m) => ({
        uniprot_id: m['Uniprot Id'],
        predicted_operator: m['Predicted operator'],
        align_score: m['Align score'],
      }))
      best.frequency_matrix = generateFrequencyMatrix(metrics)
    }
  }

  return best.consensus_score > 0 ? best : null
}
