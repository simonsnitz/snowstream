// SSE wrapper for /api/predict and helpers for /api/cache/{key}.

const decoder = new TextDecoder()

// runPipeline: kicks off /api/predict and yields parsed events.
// Usage:
//   for await (const ev of runPipeline(req, { signal })) {
//     ev.type === 'cached' | 'stage_started' | 'stage_done' | 'progress' | 'complete' | 'error'
//   }
export async function* runPipeline(req, { signal } = {}) {
  const res = await fetch('/api/predict', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(req),
    signal,
  })
  if (!res.ok) {
    throw new Error(`predict failed: HTTP ${res.status}`)
  }
  const reader = res.body.getReader()
  let buffer = ''
  while (true) {
    const { done, value } = await reader.read()
    if (done) break
    buffer += decoder.decode(value, { stream: true })
    let idx
    while ((idx = buffer.indexOf('\n\n')) !== -1) {
      const chunk = buffer.slice(0, idx)
      buffer = buffer.slice(idx + 2)
      const line = chunk.split('\n').find((l) => l.startsWith('data: '))
      if (!line) continue
      try {
        yield JSON.parse(line.slice(6))
      } catch {
        // ignore malformed
      }
    }
  }
}

export async function fetchCachedResult(key) {
  const res = await fetch(`/api/cache/${encodeURIComponent(key)}`)
  if (res.status === 404) return null
  if (!res.ok) throw new Error(`cache fetch failed: ${res.status}`)
  return await res.json()
}

export async function deleteCachedResult(key) {
  const res = await fetch(`/api/cache/${encodeURIComponent(key)}`, { method: 'DELETE' })
  if (!res.ok && res.status !== 404) throw new Error(`cache delete failed: ${res.status}`)
}
