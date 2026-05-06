// User supplies a single candidate sequence to align against all promoters.
export const params = { seq_to_align: '' }
export const label = 'Align input sequence'

export function run(_promoters, p) {
  if (!p.seq_to_align || p.seq_to_align.length < 10) return []
  return [{ seq: p.seq_to_align.toUpperCase() }]
}
