import { Alert, Button, Stack, Typography } from '@mui/material'

export default function CacheBanner({ cacheKey, cachedAt, onRerun }) {
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
