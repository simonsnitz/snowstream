// Run with: node frontend/src/operator/__tests__/verify.mjs
//
// Compares the JS operator-extraction port against fixtures dumped from the
// Python implementation. Exits non-zero on mismatch.

import assert from 'node:assert/strict'
import { readFileSync } from 'node:fs'
import { fileURLToPath } from 'node:url'
import { dirname, resolve } from 'node:path'

import { localAlign } from '../smithWaterman.js'
import {
  findOperatorInIntergenic,
  getConsensus,
  getConsensusScore,
  generateFrequencyMatrix,
} from '../consensus.js'
import { run as invertedRepeatsRun } from '../methods/invertedRepeats.js'

const here = dirname(fileURLToPath(import.meta.url))
const fixtures = JSON.parse(readFileSync(resolve(here, 'fixtures.json'), 'utf8'))

let passed = 0
let failed = 0
function test(name, fn) {
  try {
    fn()
    console.log(`  ok  ${name}`)
    passed++
  } catch (e) {
    console.error(`  FAIL ${name}`)
    console.error(`       ${e.message}`)
    failed++
  }
}

console.log('Smith-Waterman:')
for (const c of fixtures.smith_waterman) {
  test(c.case, () => {
    const got = localAlign(c.a, c.b, {
      match: c.params.match,
      mismatch: c.params.mismatch,
      gapOpen: c.params.gap_open,
      gapExtend: c.params.gap_extend,
    })
    // Bio.pairwise2 reports start/end as column indices in the padded
    // alignment strings, while our JS port reports indices in the raw `a`.
    // They diverge when there are gaps; the score is what matters for
    // downstream operator extraction. When there are tied optimal local
    // alignments, the two implementations may choose different ones.
    assert.equal(got.score, c.expected.score, 'score')
  })
}

console.log('\nfindOperatorInIntergenic:')
{
  const f = fixtures.find_operator_in_intergenic
  test('basic', () => {
    const got = findOperatorInIntergenic(f.intergenic, f.operator, f.params)
    if (f.expected === null) {
      assert.equal(got, null)
    } else {
      assert.ok(got, 'expected non-null result')
      assert.equal(got.score, f.expected.score, 'score')
      assert.equal(got.operator, f.expected.operator, 'operator string')
    }
  })
}

console.log('\nfindBestPalindrome:')
{
  const f = fixtures.find_best_palindrome
  test('produces same best operators', () => {
    const got = invertedRepeatsRun([f.intergenic], {
      min_operator_length: f.shortest,
      max_operator_length: f.longest,
      win_score: f.win_score,
      loss_score: f.loss_score,
      spacer_penalty: f.spacer_penalty,
    })
    const expectedSeqs = (f.expected || []).map((o) => o.seq).sort()
    const gotSeqs = got.map((o) => o.seq).sort()
    assert.deepEqual(gotSeqs, expectedSeqs)
    if (f.expected && f.expected.length > 0) {
      assert.equal(got[0].score, f.expected[0].score, 'top score')
    }
  })
}

console.log('\nconsensus:')
{
  const f = fixtures.get_consensus
  test('getConsensus motif_data + num_seqs', () => {
    const got = getConsensus(f.metrics)
    assert.equal(got.num_seqs, f.expected.num_seqs)
    assert.deepEqual(got.motif_data, f.expected.motif_data)
  })
  test('getConsensusScore', () => {
    const consensus = getConsensus(f.metrics)
    const got = getConsensusScore(f.score_input, consensus, 0)
    assert.equal(got, f.expected_score)
  })
  test('generateFrequencyMatrix', () => {
    const got = generateFrequencyMatrix(f.metrics)
    // Match within rounding tolerance — Python's epsilon-add rounding can differ
    // from JS by 0.01 at the last decimal.
    assert.equal(got.length, f.expected_frequency_matrix.length)
    for (let i = 0; i < got.length; i++) {
      for (let j = 0; j < 4; j++) {
        const diff = Math.abs(got[i][j] - f.expected_frequency_matrix[i][j])
        assert.ok(diff <= 0.01, `row ${i} col ${j}: got ${got[i][j]} expected ${f.expected_frequency_matrix[i][j]}`)
      }
    }
  })
}

console.log(`\n${passed} passed, ${failed} failed`)
process.exit(failed === 0 ? 0 : 1)
