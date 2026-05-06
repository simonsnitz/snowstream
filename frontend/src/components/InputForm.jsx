import {
  Box,
  FormControl,
  FormControlLabel,
  Radio,
  RadioGroup,
  TextField,
} from '@mui/material'

const PLACEHOLDERS = {
  RefSeq: 'WP_013083972.1',
  Uniprot: 'P43506',
  'Protein sequence': 'Paste a protein sequence (≥50 amino acids)',
}

export default function InputForm({ inputMethod, inputValue, onChange, disabled }) {
  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2 }}>
      <FormControl>
        <RadioGroup
          row
          value={inputMethod}
          onChange={(e) => onChange({ inputMethod: e.target.value, inputValue: '' })}
        >
          <FormControlLabel value="RefSeq" control={<Radio />} label="RefSeq" disabled={disabled} />
          <FormControlLabel value="Uniprot" control={<Radio />} label="UniProt" disabled={disabled} />
          <FormControlLabel
            value="Protein sequence"
            control={<Radio />}
            label="Protein sequence"
            disabled={disabled}
          />
        </RadioGroup>
      </FormControl>
      <TextField
        label={inputMethod}
        placeholder={PLACEHOLDERS[inputMethod]}
        value={inputValue}
        onChange={(e) => onChange({ inputValue: e.target.value })}
        multiline={inputMethod === 'Protein sequence'}
        minRows={inputMethod === 'Protein sequence' ? 4 : 1}
        fullWidth
        disabled={disabled}
      />
    </Box>
  )
}
