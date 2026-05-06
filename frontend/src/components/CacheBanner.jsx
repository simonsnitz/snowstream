import { Alert, Button, Link, Stack, Typography } from '@mui/material'

function pct(v) {
  if (v === null || v === undefined) return '?'
  return `${Math.round(v * 100)}%`
}

function characterisationLine(record) {
  const evidence = record?.evidence
  if (!evidence) return null
  if (evidence.source === 'groovdb' && evidence.url) {
    return (
      <>
        Characterized in{' '}
        <Link href={evidence.url} target="_blank" rel="noreferrer">
          groovDB ↗
        </Link>
        .
      </>
    )
  }
  if (evidence.source === 'paperblast' && evidence.papers?.length) {
    const summary = evidence.papers
      .slice(0, 3)
      .map((p) => (p.year ? `${p.year} (${p.doi})` : p.doi))
      .join(', ')
    const extra = evidence.papers.length > 3 ? ` +${evidence.papers.length - 3} more` : ''
    return <>Literature: {summary}{extra}</>
  }
  return null
}

export default function CacheBanner({ cacheKey, cachedAt, matchedVia, record, onRerun }) {
  if (matchedVia) {
    const characterisation = characterisationLine(record)
    return (
      <Alert
        severity="info"
        action={
          <Button color="inherit" size="small" onClick={onRerun}>
            Re-run with full BLAST
          </Button>
        }
      >
        <Stack spacing={0.5}>
          <Typography variant="body2">
            Matched cached representative <strong>{matchedVia.uniprot_id}</strong>{' '}
            ({pct(matchedVia.identity)} identity, {pct(matchedVia.coverage)} coverage)
            {' — '}
            family <strong>{matchedVia.family}</strong>
          </Typography>
          {characterisation && (
            <Typography variant="body2" color="text.secondary">
              {characterisation}
            </Typography>
          )}
        </Stack>
      </Alert>
    )
  }

  if (!cacheKey || !cachedAt) return null
  const when = new Date(cachedAt).toLocaleString()
  return (
    <Alert
      severity="info"
      action={
        <Button color="inherit" size="small" onClick={onRerun}>
          Re-run
        </Button>
      }
    >
      <Stack direction="row" spacing={1} alignItems="center">
        <Typography variant="body2">Showing cached result from {when}</Typography>
      </Stack>
    </Alert>
  )
}
