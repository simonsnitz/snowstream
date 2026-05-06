import { useState } from 'react'
import { Box, IconButton, Paper, Stack, Tooltip, Typography } from '@mui/material'
import ContentCopyIcon from '@mui/icons-material/ContentCopy'
import CheckIcon from '@mui/icons-material/Check'
import MotifLogo from './MotifLogo.jsx'

function CopyButton({ text }) {
  const [copied, setCopied] = useState(false)
  async function copy() {
    try {
      await navigator.clipboard.writeText(text)
      setCopied(true)
      setTimeout(() => setCopied(false), 1500)
    } catch {
      // ignore
    }
  }
  return (
    <Tooltip title={copied ? 'Copied' : 'Copy'}>
      <IconButton size="small" onClick={copy}>
        {copied ? <CheckIcon fontSize="small" /> : <ContentCopyIcon fontSize="small" />}
      </IconButton>
    </Tooltip>
  )
}

function NativePromoter({ result, homologs }) {
  const native = result.native_operator
  if (!native) return null
  const promoter = homologs.find((h) => h.promoter)?.promoter
  if (!promoter) return null
  const upper = promoter.toUpperCase()
  const idx = upper.indexOf(native.toUpperCase())
  if (idx === -1) {
    return (
      <Typography sx={{ fontFamily: 'monospace', wordBreak: 'break-all' }}>{promoter}</Typography>
    )
  }
  const before = promoter.slice(0, idx)
  const middle = promoter.slice(idx, idx + native.length)
  const after = promoter.slice(idx + native.length)
  return (
    <Typography sx={{ fontFamily: 'monospace', wordBreak: 'break-all' }}>
      <span>{before}</span>
      <span style={{ color: '#d62728', fontWeight: 700 }}>{middle}</span>
      <span>{after}</span>
    </Typography>
  )
}

export default function ResultsPanel({ result, homologs }) {
  if (!result) return null

  return (
    <Paper sx={{ p: 3 }} variant="outlined">
      <Typography variant="h5" gutterBottom>
        Results
      </Typography>
      <Stack direction={{ xs: 'column', md: 'row' }} spacing={4}>
        <Stack spacing={2}>
          <Box>
            <Typography variant="caption" color="text.secondary">
              Conservation score
            </Typography>
            <Typography variant="h4">{result.consensus_score}</Typography>
          </Box>
          <Box>
            <Typography variant="caption" color="text.secondary">
              Sequences aligned
            </Typography>
            <Typography variant="h4">{result.num_seqs}</Typography>
          </Box>
        </Stack>
        <Box sx={{ flex: 1 }}>
          <Typography variant="caption" color="text.secondary">
            Predicted promoter region
          </Typography>
          <NativePromoter result={result} homologs={homologs} />
        </Box>
      </Stack>

      <Box sx={{ mt: 4 }}>
        <Stack direction="row" alignItems="center" spacing={1}>
          <Typography variant="caption" color="text.secondary">
            Consensus sequence
          </Typography>
          <CopyButton text={result.consensus_seq} />
        </Stack>
        <Typography sx={{ fontFamily: 'monospace', fontSize: 28, letterSpacing: 1 }}>
          {result.consensus_seq}
        </Typography>
      </Box>

      <Box sx={{ mt: 3 }}>
        <Typography variant="caption" color="text.secondary">
          Conservation motif
        </Typography>
        <MotifLogo ppm={result.frequency_matrix} />
      </Box>
    </Paper>
  )
}
