// Registry of operator-extraction methods. To add a new method, drop a file
// in ./methods/ that exports { label, params, run(promoters, params) }, and
// register it here.

import * as invertedRepeats from './methods/invertedRepeats.js'
import * as alignSequence from './methods/alignSequence.js'
import * as biomsa from './methods/biomsa.js'

export const METHODS = {
  invertedRepeats,
  alignSequence,
  biomsa,
}

export const DEFAULT_METHOD = 'invertedRepeats'

export { extractOperators, DEFAULT_PARAMS } from './pipeline.js'
