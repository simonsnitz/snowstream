import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Box,
  Checkbox,
  FormControl,
  FormControlLabel,
  FormLabel,
  Radio,
  RadioGroup,
  Stack,
  TextField,
  Typography,
} from '@mui/material'
import ExpandMoreIcon from '@mui/icons-material/ExpandMore'

function NumField({ label, value, onChange, min, max, step = 1, ...rest }) {
  return (
    <TextField
      label={label}
      type="number"
      size="small"
      value={value}
      onChange={(e) => onChange(Number(e.target.value))}
      inputProps={{ min, max, step }}
      sx={{ width: 180 }}
      {...rest}
    />
  )
}

export default function AdvancedOptions({
  blast,
  promoter,
  coordinatesMethod,
  alignment,
  onChange,
}) {
  const set = (key, val) => onChange({ [key]: val })

  return (
    <Box sx={{ width: 360 }}>
      <Typography variant="h6" sx={{ px: 2, pt: 2 }}>
        Advanced options
      </Typography>

      <Accordion defaultExpanded disableGutters>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography>BLAST</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <Stack spacing={2}>
            <NumField
              label="Identity cutoff (%)"
              value={blast.ident_cutoff}
              onChange={(v) => set('blast', { ...blast, ident_cutoff: v })}
              min={30}
              max={90}
            />
            <NumField
              label="Coverage cutoff (%)"
              value={blast.cov_cutoff}
              onChange={(v) => set('blast', { ...blast, cov_cutoff: v })}
              min={50}
              max={100}
            />
            <NumField
              label="Max homologs"
              value={blast.max_homologs}
              onChange={(v) => set('blast', { ...blast, max_homologs: v })}
              min={10}
              max={100}
            />
            <FormControlLabel
              control={
                <Checkbox
                  checked={blast.filter_redundant}
                  onChange={(e) => set('blast', { ...blast, filter_redundant: e.target.checked })}
                />
              }
              label="Filter redundant"
            />
          </Stack>
        </AccordionDetails>
      </Accordion>

      <Accordion defaultExpanded disableGutters>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography>Promoter</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <Stack spacing={2}>
            <NumField
              label="Min length"
              value={promoter.min_length}
              onChange={(v) => set('promoter', { ...promoter, min_length: v })}
              min={1}
              max={500}
            />
            <NumField
              label="Max length"
              value={promoter.max_length}
              onChange={(v) => set('promoter', { ...promoter, max_length: v })}
              min={20}
              max={9000}
            />
            <FormControl>
              <FormLabel>Genome coordinates</FormLabel>
              <RadioGroup
                row
                value={coordinatesMethod}
                onChange={(e) => set('coordinatesMethod', e.target.value)}
              >
                <FormControlLabel value="batch" control={<Radio />} label="Batch" />
                <FormControlLabel value="individually" control={<Radio />} label="Per homolog" />
              </RadioGroup>
            </FormControl>
          </Stack>
        </AccordionDetails>
      </Accordion>

      <Accordion disableGutters>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          <Typography>Alignment</Typography>
        </AccordionSummary>
        <AccordionDetails>
          <Stack spacing={2}>
            <NumField
              label="Extension length"
              value={alignment.extension_length}
              onChange={(v) => set('alignment', { ...alignment, extension_length: v })}
              min={0}
              max={10}
            />
            <NumField
              label="Gap open"
              value={alignment.gap_open}
              onChange={(v) => set('alignment', { ...alignment, gap_open: v })}
              min={-999}
              max={0}
            />
            <NumField
              label="Gap extend"
              value={alignment.gap_extend}
              onChange={(v) => set('alignment', { ...alignment, gap_extend: v })}
              min={-999}
              max={0}
            />
            <NumField
              label="Match"
              value={alignment.align_match}
              onChange={(v) => set('alignment', { ...alignment, align_match: v })}
              min={1}
              max={100}
            />
            <NumField
              label="Mismatch"
              value={alignment.align_mismatch}
              onChange={(v) => set('alignment', { ...alignment, align_mismatch: v })}
              min={-100}
              max={10}
              step={0.5}
            />
          </Stack>
        </AccordionDetails>
      </Accordion>
    </Box>
  )
}
