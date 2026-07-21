// Browse view for the high-confidence curated operator dataset (4,681 TetR
// centroids where at least one algorithm produced a strong consensus).
//
// The dataset is bundled with the frontend build via
// `src/data/curated_operators.json` — no backend calls; everything is
// client-side filtering, sorting, pagination.
//
// Row click opens a detail dialog (CuratedEntryDetail.jsx) with the full
// protein sequence and side-by-side V0 vs V2.6 predictions.

import { useMemo, useState } from 'react'
import {
  Box,
  InputAdornment,
  Paper,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TablePagination,
  TableRow,
  TableSortLabel,
  TextField,
  Tooltip,
  Typography,
} from '@mui/material'
import SearchIcon from '@mui/icons-material/Search'

import curatedData from '../data/curated_operators.json'
import CuratedEntryDetail from './CuratedEntryDetail.jsx'

// One column definition = { id, label, align, getValue }.
// `getValue` is what we sort and search by. What's displayed can differ
// (e.g. we display a motif with mixed case for palindrome context, but
// sort/search on the plain uppercase string).
const COLUMNS = [
  { id: 'uniprot_id', label: 'UniProt', align: 'left',
    getValue: (r) => r.uniprot_id || '' },
  { id: 'ncbi_accession', label: 'NCBI', align: 'left',
    getValue: (r) => r.ncbi_accession || '' },
  { id: 'n_homologs', label: '# homologs', align: 'right',
    getValue: (r) => r.n_homologs ?? 0 },
  { id: 'v0_motif', label: 'V0 motif', align: 'left',
    getValue: (r) => (r.operator_v0?.motif || '').toUpperCase() },
  { id: 'v0_score', label: 'V0 cons %', align: 'right',
    getValue: (r) => r.operator_v0?.consensus_score ?? 0 },
  { id: 'v26_motif', label: 'V2.6 motif', align: 'left',
    getValue: (r) => (r.operator_v26?.motif || '').toUpperCase() },
  { id: 'v26_score', label: 'V2.6 cons %', align: 'right',
    getValue: (r) => r.operator_v26?.consensus_score ?? 0 },
  { id: 'v26_rerank', label: 'V2.6 rerank', align: 'right',
    getValue: (r) => r.operator_v26?.rerank_score ?? 0 },
]

// Motifs come with lowercase flanking + uppercase palindrome core. Render
// with a very light uppercase / bold uppercase mix so the operator "arms"
// pop visually. Keep the whole thing monospace.
function MotifCell({ motif }) {
  if (!motif) return <em>—</em>
  return (
    <Box component="span" sx={{
      fontFamily: 'ui-monospace, "SF Mono", Menlo, Consolas, monospace',
      fontSize: 13,
      whiteSpace: 'nowrap',
    }}>
      {motif.split('').map((ch, i) => (
        <span
          key={i}
          style={{
            color: ch === ch.toUpperCase() && ch !== '-' ? '#111' : '#999',
            fontWeight: ch === ch.toUpperCase() && ch !== '-' ? 600 : 400,
          }}
        >
          {ch}
        </span>
      ))}
    </Box>
  )
}

function fmtScore(v) {
  if (v == null) return '—'
  const n = Number(v)
  if (!Number.isFinite(n)) return '—'
  return n.toFixed(1)
}

function fmtRerank(v) {
  if (v == null) return '—'
  const n = Number(v)
  if (!Number.isFinite(n)) return '—'
  return n.toFixed(3)
}

export default function CuratedSet() {
  const [search, setSearch] = useState('')
  const [orderBy, setOrderBy] = useState('v26_score')
  const [order, setOrder] = useState('desc')
  const [page, setPage] = useState(0)
  const [rowsPerPage, setRowsPerPage] = useState(25)
  const [openEntry, setOpenEntry] = useState(null)

  // Precompute lowercase search haystacks once per row so the filter is fast
  // even at 4,681 rows.
  const rows = useMemo(() => {
    return curatedData.map((r) => ({
      ...r,
      _search: [
        r.uniprot_id,
        r.ncbi_accession,
        (r.operator_v0?.motif || ''),
        (r.operator_v26?.motif || ''),
        r.genome,
      ].filter(Boolean).join(' ').toLowerCase(),
    }))
  }, [])

  const filtered = useMemo(() => {
    const q = search.trim().toLowerCase()
    if (!q) return rows
    return rows.filter((r) => r._search.includes(q))
  }, [rows, search])

  const sorted = useMemo(() => {
    const col = COLUMNS.find((c) => c.id === orderBy) || COLUMNS[0]
    const cmp = (a, b) => {
      const av = col.getValue(a)
      const bv = col.getValue(b)
      if (av < bv) return order === 'asc' ? -1 : 1
      if (av > bv) return order === 'asc' ? 1 : -1
      return 0
    }
    // Copy so we don't mutate the memoised filtered array.
    return [...filtered].sort(cmp)
  }, [filtered, orderBy, order])

  const paged = useMemo(() => {
    const start = page * rowsPerPage
    return sorted.slice(start, start + rowsPerPage)
  }, [sorted, page, rowsPerPage])

  function handleSort(colId) {
    if (orderBy === colId) {
      setOrder((prev) => (prev === 'asc' ? 'desc' : 'asc'))
    } else {
      setOrderBy(colId)
      // Sensible default: numeric columns descend on first click (highest
      // scores at the top); text columns ascend.
      const col = COLUMNS.find((c) => c.id === colId)
      setOrder(typeof col?.getValue({}) === 'number' ? 'desc' : 'asc')
    }
    setPage(0)
  }

  return (
    <Stack spacing={2}>
      <Paper sx={{ p: 2 }} variant="outlined">
        <Stack direction={{ xs: 'column', md: 'row' }} spacing={2}
                alignItems={{ md: 'center' }} justifyContent="space-between">
          <Box>
            <Typography variant="h6">TetR high-confidence set</Typography>
            <Typography variant="body2" color="text.secondary">
              {curatedData.length.toLocaleString()} cluster centroids with
              ≥ 5 homologs and consensus &gt; 75 % in V0 and/or V2.6.
              Click a row to see the protein sequence and both algorithm
              predictions in detail.
            </Typography>
          </Box>
          <TextField
            size="small"
            placeholder="Search UniProt, NCBI, motif, genome…"
            value={search}
            onChange={(e) => { setSearch(e.target.value); setPage(0) }}
            sx={{ minWidth: 300 }}
            InputProps={{
              startAdornment: (
                <InputAdornment position="start">
                  <SearchIcon fontSize="small" />
                </InputAdornment>
              ),
            }}
          />
        </Stack>
      </Paper>

      <Paper variant="outlined">
        <TableContainer sx={{ maxHeight: '70vh' }}>
          <Table size="small" stickyHeader>
            <TableHead>
              <TableRow>
                {COLUMNS.map((col) => (
                  <TableCell
                    key={col.id}
                    align={col.align}
                    sortDirection={orderBy === col.id ? order : false}
                    sx={{ fontWeight: 700, whiteSpace: 'nowrap' }}
                  >
                    <TableSortLabel
                      active={orderBy === col.id}
                      direction={orderBy === col.id ? order : 'asc'}
                      onClick={() => handleSort(col.id)}
                    >
                      {col.label}
                    </TableSortLabel>
                  </TableCell>
                ))}
              </TableRow>
            </TableHead>
            <TableBody>
              {paged.map((r) => (
                <TableRow
                  key={r.uniprot_id}
                  hover
                  onClick={() => setOpenEntry(r)}
                  sx={{ cursor: 'pointer' }}
                >
                  <TableCell sx={{ fontFamily: 'ui-monospace, monospace', fontSize: 13 }}>
                    {r.uniprot_id}
                  </TableCell>
                  <TableCell sx={{ fontFamily: 'ui-monospace, monospace', fontSize: 13 }}>
                    {r.ncbi_accession || <em style={{ color: '#999' }}>—</em>}
                  </TableCell>
                  <TableCell align="right">{r.n_homologs}</TableCell>
                  <TableCell>
                    <Tooltip title={r.operator_v0?.motif || ''} enterDelay={400}>
                      <span><MotifCell motif={r.operator_v0?.motif} /></span>
                    </Tooltip>
                  </TableCell>
                  <TableCell align="right">{fmtScore(r.operator_v0?.consensus_score)}</TableCell>
                  <TableCell>
                    <Tooltip title={r.operator_v26?.motif || ''} enterDelay={400}>
                      <span><MotifCell motif={r.operator_v26?.motif} /></span>
                    </Tooltip>
                  </TableCell>
                  <TableCell align="right">{fmtScore(r.operator_v26?.consensus_score)}</TableCell>
                  <TableCell align="right">{fmtRerank(r.operator_v26?.rerank_score)}</TableCell>
                </TableRow>
              ))}
              {paged.length === 0 && (
                <TableRow>
                  <TableCell colSpan={COLUMNS.length} align="center">
                    <Typography variant="body2" color="text.secondary" sx={{ py: 3 }}>
                      No entries match "{search}".
                    </Typography>
                  </TableCell>
                </TableRow>
              )}
            </TableBody>
          </Table>
        </TableContainer>
        <TablePagination
          rowsPerPageOptions={[25, 50, 100, 200]}
          component="div"
          count={sorted.length}
          rowsPerPage={rowsPerPage}
          page={page}
          onPageChange={(_e, p) => setPage(p)}
          onRowsPerPageChange={(e) => {
            setRowsPerPage(parseInt(e.target.value, 10))
            setPage(0)
          }}
        />
      </Paper>

      <CuratedEntryDetail
        entry={openEntry}
        onClose={() => setOpenEntry(null)}
      />
    </Stack>
  )
}
