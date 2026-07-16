// Detail dialog for one row of the CuratedSet browse table.
//
// Shows: identifiers, genome, homolog count, protein sequence (with copy
// button), and the two algorithm predictions side-by-side.

import { useState } from 'react'
import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Alert,
  Box,
  Button,
  Chip,
  Dialog,
  DialogContent,
  DialogTitle,
  Divider,
  IconButton,
  Link,
  Paper,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Typography,
} from '@mui/material'
import CloseIcon from '@mui/icons-material/Close'
import ContentCopyIcon from '@mui/icons-material/ContentCopy'
import ExpandMoreIcon from '@mui/icons-material/ExpandMore'
import OpenInNewIcon from '@mui/icons-material/OpenInNew'

import MotifLogo from './MotifLogo.jsx'

function MotifDisplay({ motif }) {
  if (!motif) return <em style={{ color: '#999' }}>none</em>
  return (
    <Box component="span" sx={{
      fontFamily: 'ui-monospace, "SF Mono", Menlo, Consolas, monospace',
      fontSize: 15,
      letterSpacing: '0.5px',
      wordBreak: 'break-all',
    }}>
      {motif.split('').map((ch, i) => {
        const isCore = ch === ch.toUpperCase() && ch !== '-'
        return (
          <span
            key={i}
            style={{
              color: isCore ? '#111' : '#999',
              fontWeight: isCore ? 700 : 400,
            }}
          >
            {ch}
          </span>
        )
      })}
    </Box>
  )
}

function VersionPanel({ label, block, extra }) {
  if (!block) return null
  const {
    motif,
    consensus_score: cons,
    rerank_score: rerank,
    native_for_centroid: native,
    frequency_matrix: ppm,
    aligned_operators: aligned,
  } = block
  return (
    <Paper variant="outlined" sx={{ p: 2, flex: 1 }}>
      <Stack spacing={1.5}>
        <Stack direction="row" spacing={1} alignItems="center">
          <Typography variant="h6" sx={{ fontWeight: 700 }}>{label}</Typography>
          {extra}
        </Stack>

        <Box>
          <Typography variant="caption" color="text.secondary">
            Consensus motif
          </Typography>
          <Box sx={{ mt: 0.5 }}>
            <MotifDisplay motif={motif} />
          </Box>
        </Box>

        {ppm?.length > 0 && (
          <Box>
            <Typography variant="caption" color="text.secondary">
              Motif logo (PPM from {aligned?.length ?? 0} aligned operators)
            </Typography>
            <Box sx={{ mt: 0.5, minHeight: 80 }}>
              <MotifLogo ppm={ppm} />
            </Box>
          </Box>
        )}

        <Stack direction="row" spacing={3}>
          <Box>
            <Typography variant="caption" color="text.secondary">
              Consensus score
            </Typography>
            <Typography variant="h6" sx={{ fontWeight: 600 }}>
              {cons != null ? Number(cons).toFixed(1) : '—'}
            </Typography>
          </Box>
          {rerank != null && (
            <Box>
              <Typography variant="caption" color="text.secondary">
                Rerank score
              </Typography>
              <Typography variant="h6" sx={{ fontWeight: 600 }}>
                {Number(rerank).toFixed(3)}
              </Typography>
            </Box>
          )}
        </Stack>

        {native && (
          <Box>
            <Typography variant="caption" color="text.secondary">
              Operator in centroid&apos;s own promoter
            </Typography>
            <Box sx={{ mt: 0.5 }}>
              <MotifDisplay motif={native} />
            </Box>
          </Box>
        )}

        {aligned?.length > 0 && <AlignedOperatorsAccordion aligned={aligned} />}
      </Stack>
    </Paper>
  )
}


// Expandable table of the per-homolog aligned operators that fed the
// consensus. Sorted by align_score desc so the strongest hits are visible
// without scrolling. Capped visual height with an inner scroll — the list
// can be up to 100 entries per version.
function AlignedOperatorsAccordion({ aligned }) {
  const rows = [...aligned].sort((a, b) => {
    const as = a.align_score ?? -Infinity
    const bs = b.align_score ?? -Infinity
    return bs - as
  })
  return (
    <Accordion
      disableGutters
      elevation={0}
      sx={{
        border: '1px solid',
        borderColor: 'divider',
        '&:before': { display: 'none' },
      }}
    >
      <AccordionSummary expandIcon={<ExpandMoreIcon />}>
        <Typography variant="body2">
          Per-homolog aligned operators ({rows.length})
        </Typography>
      </AccordionSummary>
      <AccordionDetails sx={{ p: 0 }}>
        <TableContainer sx={{ maxHeight: 340 }}>
          <Table size="small" stickyHeader>
            <TableHead>
              <TableRow>
                <TableCell sx={{ fontWeight: 700 }}>UniProt</TableCell>
                <TableCell sx={{ fontWeight: 700 }}>Predicted operator</TableCell>
                <TableCell align="right" sx={{ fontWeight: 700 }}>Align score</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {rows.map((r, i) => (
                <TableRow key={`${r.uniprot_id}-${i}`} hover>
                  <TableCell sx={{ fontFamily: 'ui-monospace, monospace', fontSize: 12, whiteSpace: 'nowrap' }}>
                    {r.uniprot_id ? (
                      <Link
                        href={`https://www.uniprot.org/uniprotkb/${r.uniprot_id}`}
                        target="_blank"
                        rel="noreferrer"
                        underline="hover"
                      >
                        {r.uniprot_id}
                      </Link>
                    ) : <em style={{ color: '#999' }}>—</em>}
                  </TableCell>
                  <TableCell><MotifDisplay motif={r.operator} /></TableCell>
                  <TableCell align="right" sx={{ fontSize: 12 }}>
                    {r.align_score != null ? Number(r.align_score).toFixed(1) : '—'}
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </TableContainer>
      </AccordionDetails>
    </Accordion>
  )
}

// Format a protein sequence into fixed-width 10-aa blocks so the alignment
// eye-check is easy. Renders "  1 MKTLA VSTKD LGSMR..." with position
// counters every 60 chars, matching the FASTA / clustal convention.
function FormattedSequence({ seq }) {
  if (!seq) return null
  const blocks = []
  const perLine = 60
  const perBlock = 10
  for (let i = 0; i < seq.length; i += perLine) {
    const line = seq.slice(i, i + perLine)
    const grouped = []
    for (let j = 0; j < line.length; j += perBlock) {
      grouped.push(line.slice(j, j + perBlock))
    }
    blocks.push(
      <Box key={i} sx={{ display: 'flex' }}>
        <Box sx={{ width: 60, color: 'text.secondary', textAlign: 'right', pr: 2 }}>
          {i + 1}
        </Box>
        <Box>{grouped.join(' ')}</Box>
      </Box>
    )
  }
  return (
    <Box sx={{
      fontFamily: 'ui-monospace, "SF Mono", Menlo, Consolas, monospace',
      fontSize: 12,
      lineHeight: 1.6,
      whiteSpace: 'pre',
    }}>
      {blocks}
    </Box>
  )
}

export default function CuratedEntryDetail({ entry, onClose }) {
  const [copied, setCopied] = useState(false)

  if (!entry) return null

  async function copySequence() {
    try {
      await navigator.clipboard.writeText(entry.protein_sequence || '')
      setCopied(true)
      setTimeout(() => setCopied(false), 1500)
    } catch {
      // clipboard may be unavailable in some contexts — silent fail
    }
  }

  const uniprotURL = `https://www.uniprot.org/uniprotkb/${entry.uniprot_id}`
  const ncbiURL = entry.ncbi_accession
    ? `https://www.ncbi.nlm.nih.gov/protein/${entry.ncbi_accession}`
    : null

  return (
    <Dialog
      open={!!entry}
      onClose={onClose}
      fullWidth
      maxWidth="lg"
      scroll="paper"
    >
      <DialogTitle sx={{ pr: 6 }}>
        <Stack direction="row" spacing={1} alignItems="baseline">
          <Typography variant="h6">{entry.uniprot_id}</Typography>
          {entry.ncbi_accession && (
            <Typography variant="body2" color="text.secondary">
              · {entry.ncbi_accession}
            </Typography>
          )}
        </Stack>
        <IconButton
          onClick={onClose}
          sx={{ position: 'absolute', right: 8, top: 8 }}
          aria-label="close"
        >
          <CloseIcon />
        </IconButton>
      </DialogTitle>

      <DialogContent dividers>
        <Stack spacing={3}>
          {/* Identifiers + genome + homolog count */}
          <Stack direction={{ xs: 'column', md: 'row' }} spacing={2}>
            <Paper variant="outlined" sx={{ p: 2, flex: 1 }}>
              <Typography variant="caption" color="text.secondary">UniProt</Typography>
              <Box sx={{ mt: 0.5 }}>
                <Link href={uniprotURL} target="_blank" rel="noreferrer" underline="hover">
                  {entry.uniprot_id}
                  <OpenInNewIcon sx={{ fontSize: 14, ml: 0.5, verticalAlign: 'middle' }} />
                </Link>
              </Box>
            </Paper>
            <Paper variant="outlined" sx={{ p: 2, flex: 1 }}>
              <Typography variant="caption" color="text.secondary">NCBI RefSeq</Typography>
              <Box sx={{ mt: 0.5 }}>
                {ncbiURL ? (
                  <Link href={ncbiURL} target="_blank" rel="noreferrer" underline="hover">
                    {entry.ncbi_accession}
                    <OpenInNewIcon sx={{ fontSize: 14, ml: 0.5, verticalAlign: 'middle' }} />
                  </Link>
                ) : (
                  <em style={{ color: '#999' }}>not linked</em>
                )}
              </Box>
            </Paper>
            <Paper variant="outlined" sx={{ p: 2, flex: 1 }}>
              <Typography variant="caption" color="text.secondary">Genome</Typography>
              <Box sx={{ mt: 0.5, fontFamily: 'ui-monospace, monospace', fontSize: 13 }}>
                {entry.genome || <em style={{ color: '#999' }}>—</em>}
              </Box>
            </Paper>
            <Paper variant="outlined" sx={{ p: 2, flex: 1 }}>
              <Typography variant="caption" color="text.secondary">Homologs used</Typography>
              <Typography variant="h6" sx={{ mt: 0.5, fontWeight: 600 }}>
                {entry.n_homologs}
              </Typography>
            </Paper>
          </Stack>

          {/* Side-by-side V0 / V2.6 predictions */}
          <Box>
            <Typography variant="overline" color="text.secondary">
              Algorithm predictions
            </Typography>
            <Stack direction={{ xs: 'column', md: 'row' }} spacing={2} sx={{ mt: 1 }}>
              <VersionPanel
                label="V0"
                block={entry.operator_v0}
                extra={<Chip size="small" label="legacy" />}
              />
              <VersionPanel
                label="V2.6"
                block={entry.operator_v26}
                extra={<Chip size="small" color="primary" label="widened + gc_medium_weak" />}
              />
            </Stack>
          </Box>

          {/* Protein sequence */}
          <Box>
            <Stack direction="row" spacing={1} alignItems="center" sx={{ mb: 1 }}>
              <Typography variant="overline" color="text.secondary" sx={{ flexGrow: 1 }}>
                Protein sequence ({(entry.protein_sequence || '').length} aa)
              </Typography>
              <Button
                size="small"
                variant="outlined"
                startIcon={<ContentCopyIcon />}
                onClick={copySequence}
              >
                {copied ? 'Copied!' : 'Copy'}
              </Button>
            </Stack>
            <Paper variant="outlined" sx={{ p: 2, bgcolor: 'grey.50', overflowX: 'auto' }}>
              {entry.protein_sequence ? (
                <FormattedSequence seq={entry.protein_sequence} />
              ) : (
                <Alert severity="warning">No sequence available in members.fasta.</Alert>
              )}
            </Paper>
          </Box>

          <Divider />

          <Typography variant="caption" color="text.secondary">
            <strong>Motif convention:</strong> lowercase = flanking bases (extension
            context, 5 bp on each side); UPPERCASE = the palindromic core the
            algorithm identified as the operator.
          </Typography>
        </Stack>
      </DialogContent>
    </Dialog>
  )
}
