import { useEffect, useState } from 'react'
import {
  Alert,
  Box,
  Button,
  Container,
  Drawer,
  IconButton,
  Paper,
  Stack,
  Toolbar,
  Typography,
} from '@mui/material'
import TuneIcon from '@mui/icons-material/Tune'

import InputForm from './components/InputForm.jsx'
import AdvancedOptions from './components/AdvancedOptions.jsx'
import OperatorMethodControls from './components/OperatorMethodControls.jsx'
import PipelinePanel from './components/PipelinePanel.jsx'
import ResultsPanel from './components/ResultsPanel.jsx'
import CacheBanner from './components/CacheBanner.jsx'
import { runPipeline, fetchCachedResult } from './api.js'
import { METHODS, DEFAULT_METHOD, extractOperators, DEFAULT_PARAMS } from './operator/index.js'

const DEFAULT_BLAST = { ident_cutoff: 40, cov_cutoff: 90, max_homologs: 30, filter_redundant: true }
const DEFAULT_PROMOTER = { min_length: 80, max_length: 800 }

const initialPipeline = () => ({
  active: false,
  stages: {},
  progress: null,
  error: null,
})

export default function App() {
  const [inputMethod, setInputMethod] = useState('RefSeq')
  const [inputValue, setInputValue] = useState('WP_013083972.1')
  const [blast, setBlast] = useState(DEFAULT_BLAST)
  const [promoter, setPromoter] = useState(DEFAULT_PROMOTER)
  const [coordinatesMethod, setCoordinatesMethod] = useState('batch')
  const [database, setDatabase] = useState('local_diamond')
  const [operatorMethod, setOperatorMethod] = useState(DEFAULT_METHOD)
  // Keep params per-method so user input (e.g. "Sequence to align") survives
  // toggling between methods.
  const [paramsByMethod, setParamsByMethod] = useState(() =>
    Object.fromEntries(
      Object.entries(METHODS).map(([k, m]) => [k, { ...m.params }]),
    ),
  )
  const operatorParams = paramsByMethod[operatorMethod]
  const setOperatorParams = (next) =>
    setParamsByMethod((prev) => ({ ...prev, [operatorMethod]: next }))
  const [alignment, setAlignment] = useState(DEFAULT_PARAMS)

  const [drawerOpen, setDrawerOpen] = useState(false)
  const [pipeline, setPipeline] = useState(initialPipeline)
  const [predictResult, setPredictResult] = useState(null) // homologs+protein_info from /api/predict
  const [cacheMeta, setCacheMeta] = useState(null) // { cache_key, cached_at, fromCache }

  // Hydrate from URL on initial mount: ?q=<cache_key>.
  useEffect(() => {
    const params = new URLSearchParams(window.location.search)
    const key = params.get('q')
    if (!key) return
    fetchCachedResult(key).then((hit) => {
      if (!hit) return
      setPredictResult(hit)
      setCacheMeta({ cache_key: hit.cache_key, cached_at: hit.cached_at, fromCache: true })
      // Repopulate input fields so the user can re-run with the same params.
      if (hit.input) {
        setInputMethod(hit.input.input_method)
        setInputValue(hit.input.input_value)
        if (hit.input.blast_params) setBlast(hit.input.blast_params)
        if (hit.input.promoter_params) setPromoter(hit.input.promoter_params)
        if (hit.input.get_coordinates_method) setCoordinatesMethod(hit.input.get_coordinates_method)
        if (hit.input.database) setDatabase(hit.input.database)
      }
    })
  }, [])

  // Sync URL ?q= when cache_key changes.
  useEffect(() => {
    const url = new URL(window.location.href)
    if (cacheMeta?.cache_key) {
      url.searchParams.set('q', cacheMeta.cache_key)
    } else {
      url.searchParams.delete('q')
    }
    window.history.replaceState({}, '', url)
  }, [cacheMeta?.cache_key])

  // Compute operator extraction whenever the predict result OR the operator
  // method/params/alignment params change. Methods may be async (BioMSA),
  // so this runs in an effect rather than useMemo.
  const [operatorResult, setOperatorResult] = useState(null)
  const [operatorComputing, setOperatorComputing] = useState(false)
  useEffect(() => {
    if (!predictResult?.homologs?.length) {
      setOperatorResult(null)
      return
    }
    let cancelled = false
    setOperatorComputing(true)
    Promise.resolve(
      METHODS[operatorMethod].run(
        predictResult.homologs.map((h) => h.promoter),
        operatorParams,
      ),
    )
      .then((candidates) => extractOperators(predictResult.homologs, candidates, alignment))
      .then((result) => {
        if (cancelled) return
        setOperatorResult(result)
        setOperatorComputing(false)
      })
      .catch((e) => {
        if (cancelled) return
        console.error('operator extraction failed:', e)
        setOperatorResult(null)
        setOperatorComputing(false)
      })
    return () => {
      cancelled = true
    }
  }, [predictResult, operatorMethod, operatorParams, alignment])

  async function handleSubmit({ force = false } = {}) {
    if (!inputValue) return
    setPipeline({ ...initialPipeline(), active: true, stages: { blast: 'running' } })
    setPredictResult(null)
    setCacheMeta(null)

    const req = {
      input_method: inputMethod,
      input_value: inputValue,
      blast_params: blast,
      promoter_params: promoter,
      get_coordinates_method: coordinatesMethod,
      database,
      force,
    }

    try {
      for await (const ev of runPipeline(req)) {
        if (ev.type === 'cached') {
          setPredictResult(ev.result)
          setCacheMeta({
            cache_key: ev.cache_key,
            cached_at: ev.result?.cached_at,
            fromCache: true,
          })
          setPipeline({ active: false, stages: {}, progress: null, error: null })
        } else if (ev.type === 'stage_started') {
          setPipeline((p) => ({
            ...p,
            stages: { ...p.stages, [ev.stage]: 'running' },
            progress: null,
          }))
        } else if (ev.type === 'progress') {
          setPipeline((p) => ({ ...p, progress: ev }))
        } else if (ev.type === 'stage_done') {
          setPipeline((p) => ({
            ...p,
            stages: { ...p.stages, [ev.stage]: 'done' },
            progress: null,
          }))
        } else if (ev.type === 'complete') {
          setPredictResult(ev.result)
          setCacheMeta({
            cache_key: ev.cache_key,
            cached_at: ev.result?.cached_at,
            fromCache: false,
          })
          setPipeline((p) => ({ ...p, active: false, progress: null }))
        } else if (ev.type === 'error') {
          setPipeline((p) => ({
            ...p,
            active: false,
            stages: { ...p.stages, [ev.stage]: 'error' },
            error: ev,
          }))
        }
      }
    } catch (e) {
      setPipeline((p) => ({ ...p, active: false, error: { message: String(e) } }))
    }
  }

  return (
    <>
      <Toolbar sx={{ borderBottom: 1, borderColor: 'divider' }}>
        <Typography variant="h6" sx={{ flexGrow: 1 }}>
          Snowprint
        </Typography>
        <IconButton onClick={() => setDrawerOpen(true)} aria-label="advanced options">
          <TuneIcon />
        </IconButton>
      </Toolbar>

      <Drawer anchor="right" open={drawerOpen} onClose={() => setDrawerOpen(false)}>
        <AdvancedOptions
          blast={blast}
          promoter={promoter}
          coordinatesMethod={coordinatesMethod}
          database={database}
          alignment={alignment}
          onClose={() => setDrawerOpen(false)}
          onChange={(patch) => {
            if (patch.blast) setBlast(patch.blast)
            if (patch.promoter) setPromoter(patch.promoter)
            if (patch.coordinatesMethod) setCoordinatesMethod(patch.coordinatesMethod)
            if (patch.database) setDatabase(patch.database)
            if (patch.alignment) setAlignment(patch.alignment)
          }}
        />
      </Drawer>

      <Container maxWidth="lg" sx={{ py: 4 }}>
        <Stack spacing={3}>
          <Typography variant="h4" align="center">
            Predict a regulator&apos;s DNA binding sequence
          </Typography>

          <Paper sx={{ p: 3 }} variant="outlined">
            <Stack spacing={2}>
              <InputForm
                inputMethod={inputMethod}
                inputValue={inputValue}
                onChange={(p) => {
                  if (p.inputMethod !== undefined) setInputMethod(p.inputMethod)
                  if (p.inputValue !== undefined) setInputValue(p.inputValue)
                }}
                disabled={pipeline.active}
              />
              <Box sx={{ display: 'flex', justifyContent: 'center' }}>
                <Button
                  variant="contained"
                  onClick={handleSubmit}
                  disabled={pipeline.active || !inputValue}
                >
                  {pipeline.active ? 'Running…' : 'Submit'}
                </Button>
              </Box>
            </Stack>
          </Paper>

          {cacheMeta?.fromCache && (
            <CacheBanner
              cacheKey={cacheMeta.cache_key}
              cachedAt={cacheMeta.cached_at}
              onRerun={() => handleSubmit({ force: true })}
            />
          )}

          {predictResult?.protein_info && (
            <Paper sx={{ p: 2 }} variant="outlined">
              <Typography variant="subtitle2" color="text.secondary">
                Input
              </Typography>
              <Typography>Annotation: {predictResult.protein_info.annotation}</Typography>
              <Typography>Organism: {predictResult.protein_info.organism}</Typography>
              <Typography variant="body2" color="text.secondary">
                Lineage: {predictResult.protein_info.lineage?.join(', ')}
              </Typography>
            </Paper>
          )}

          <PipelinePanel pipeline={pipeline} />

          {predictResult && !pipeline.active && !operatorComputing && !operatorResult && (
            <Alert severity="info">
              Pipeline complete. No operator candidates found with the selected method.
            </Alert>
          )}

          {operatorComputing && <Alert severity="info">Computing alignment…</Alert>}

          <ResultsPanel result={operatorResult} homologs={predictResult?.homologs || []} />

          {predictResult && (
            <OperatorMethodControls
              operatorMethod={operatorMethod}
              operatorParams={operatorParams}
              onMethodChange={setOperatorMethod}
              onParamsChange={setOperatorParams}
            />
          )}
        </Stack>
      </Container>
    </>
  )
}
