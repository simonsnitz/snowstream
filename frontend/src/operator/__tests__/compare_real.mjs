// Run the JS operator pipeline on a cached homolog set and print the result,
// for comparison against running Python fetch_operator on the same data.
//
// Usage: node frontend/src/operator/__tests__/compare_real.mjs cache/<key>.json
import { readFileSync } from 'node:fs'
import { extractOperators, METHODS, DEFAULT_PARAMS } from '../index.js'

const cachePath = process.argv[2]
if (!cachePath) {
  console.error('usage: node compare_real.mjs <cache-file>')
  process.exit(1)
}
const cache = JSON.parse(readFileSync(cachePath, 'utf8'))
const homologs = cache.homologs.filter((h) => h.promoter)
console.log('homologs with promoter:', homologs.length)

const candidates = METHODS.invertedRepeats.run(
  homologs.map((h) => h.promoter),
  METHODS.invertedRepeats.params,
)
console.log('candidates:', candidates.length)

const result = extractOperators(homologs, candidates, DEFAULT_PARAMS)
if (!result) {
  console.log('no result')
  process.exit(0)
}
console.log('js num_seqs:', result.num_seqs)
console.log('js consensus_score:', result.consensus_score)
console.log('js consensus_seq:', result.consensus_seq)
console.log('js native_operator:', result.native_operator)
