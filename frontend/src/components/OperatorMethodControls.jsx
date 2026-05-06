import {
  FormControl,
  FormControlLabel,
  FormLabel,
  Paper,
  Radio,
  RadioGroup,
  Stack,
  TextField,
  Typography,
} from '@mui/material'

import { METHODS } from '../operator/index.js'

function NumField({ label, value, onChange, min, max, step = 1 }) {
  return (
    <TextField
      label={label}
      type="number"
      size="small"
      value={value}
      onChange={(e) => onChange(Number(e.target.value))}
      inputProps={{ min, max, step }}
      sx={{ width: 180 }}
    />
  )
}

export default function OperatorMethodControls({
  operatorMethod,
  operatorParams,
  onMethodChange,
  onParamsChange,
}) {
  return (
    <Paper sx={{ p: 2 }} variant="outlined">
      <Typography variant="subtitle2" gutterBottom>
        Search method
      </Typography>
      <Stack
        direction={{ xs: 'column', md: 'row' }}
        spacing={3}
        alignItems={{ md: 'flex-start' }}
      >
        <FormControl>
          <FormLabel sx={{ fontSize: 12 }}>Method</FormLabel>
          <RadioGroup
            value={operatorMethod}
            onChange={(e) => onMethodChange(e.target.value)}
          >
            {Object.entries(METHODS).map(([key, m]) => (
              <FormControlLabel
                key={key}
                value={key}
                control={<Radio size="small" />}
                label={m.label}
              />
            ))}
          </RadioGroup>
        </FormControl>

        {operatorMethod === 'invertedRepeats' && (
          <Stack direction="row" spacing={2} flexWrap="wrap" useFlexGap>
            <NumField
              label="Match score"
              value={operatorParams.win_score}
              onChange={(v) => onParamsChange({ ...operatorParams, win_score: v })}
              min={0}
              max={10}
            />
            <NumField
              label="Mismatch score"
              value={operatorParams.loss_score}
              onChange={(v) => onParamsChange({ ...operatorParams, loss_score: v })}
              min={-10}
              max={0}
            />
            <NumField
              label="Min operator length"
              value={operatorParams.min_operator_length}
              onChange={(v) =>
                onParamsChange({ ...operatorParams, min_operator_length: v })
              }
              min={3}
              max={10}
            />
            <NumField
              label="Max operator length"
              value={operatorParams.max_operator_length}
              onChange={(v) =>
                onParamsChange({ ...operatorParams, max_operator_length: v })
              }
              min={11}
              max={40}
            />
          </Stack>
        )}

        {operatorMethod === 'biomsa' && (
          <Stack direction="row" spacing={2} flexWrap="wrap" useFlexGap>
            <NumField
              label="Min operator length"
              value={operatorParams.min_operator_length}
              onChange={(v) =>
                onParamsChange({ ...operatorParams, min_operator_length: v })
              }
              min={5}
              max={40}
            />
            <NumField
              label="Max operator length"
              value={operatorParams.max_operator_length}
              onChange={(v) =>
                onParamsChange({ ...operatorParams, max_operator_length: v })
              }
              min={5}
              max={60}
            />
          </Stack>
        )}

        {operatorMethod === 'alignSequence' && (
          <TextField
            label="Sequence to align"
            size="small"
            value={operatorParams.seq_to_align || ''}
            onChange={(e) =>
              onParamsChange({ ...operatorParams, seq_to_align: e.target.value })
            }
            sx={{ minWidth: 320 }}
          />
        )}
      </Stack>
    </Paper>
  )
}
