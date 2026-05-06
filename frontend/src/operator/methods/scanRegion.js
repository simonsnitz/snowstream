// Slide a 24bp window across the reference promoter, step 12.
export const params = {}
export const label = 'Scan promoter region'

export function run(promoters) {
  const reference = promoters.find((s) => s != null)
  if (!reference) return []
  const out = []
  for (let c = 0; c < reference.length - 24; c += 12) {
    out.push({ seq: reference.slice(c, c + 24) })
  }
  return out
}
