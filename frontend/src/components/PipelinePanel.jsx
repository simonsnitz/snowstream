import { Alert, Box, LinearProgress, Paper, Stack, Typography } from '@mui/material'

const STAGES = [
  { key: 'blast', label: 'BLAST' },
  { key: 'coordinates', label: 'Genome coordinates' },
  { key: 'operons', label: 'Operons + promoters' },
]

function StageRow({ label, status, progress }) {
  const colorByStatus = {
    pending: 'text.disabled',
    running: 'primary.main',
    done: 'success.main',
    error: 'error.main',
  }
  return (
    <Box>
      <Typography sx={{ color: colorByStatus[status] }}>
        {label}
        {status === 'running' && progress?.total
          ? ` — ${progress.current}/${progress.total}${progress.uniprot_id ? ` (${progress.uniprot_id})` : ''}`
          : ''}
        {status === 'running' && progress?.elapsed_seconds !== undefined
          ? ` — ${progress.message || 'in progress'} (${progress.elapsed_seconds}s elapsed)`
          : ''}
        {status === 'done' && ' ✓'}
      </Typography>
      {status === 'running' && (
        <LinearProgress
          variant={progress?.total ? 'determinate' : 'indeterminate'}
          value={progress?.total ? Math.round((progress.current / Math.max(progress.total, 1)) * 100) : 0}
        />
      )}
    </Box>
  )
}

export default function PipelinePanel({ pipeline }) {
  if (!pipeline.active && !pipeline.error) return null

  return (
    <Paper sx={{ p: 2 }} variant="outlined">
      <Stack spacing={1.5}>
        {STAGES.map((s) => (
          <StageRow
            key={s.key}
            label={s.label}
            status={pipeline.stages[s.key] || 'pending'}
            progress={pipeline.progress?.stage === s.key ? pipeline.progress : null}
          />
        ))}
        {pipeline.error && <Alert severity="error">{pipeline.error.message}</Alert>}
      </Stack>
    </Paper>
  )
}
