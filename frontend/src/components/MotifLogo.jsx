import { useRef } from 'react'
import { Box, Button, Stack } from '@mui/material'
import DownloadIcon from '@mui/icons-material/Download'
import { DNALogo } from 'logojs-react'

export default function MotifLogo({ ppm }) {
  const containerRef = useRef(null)
  if (!ppm?.length) return null

  function downloadPNG() {
    const svg = containerRef.current?.querySelector('svg')
    if (!svg) return
    const xml = new XMLSerializer().serializeToString(svg)
    const blob = new Blob([xml], { type: 'image/svg+xml;charset=utf-8' })
    const url = URL.createObjectURL(blob)
    const img = new Image()
    img.onload = () => {
      const scale = 3
      const rect = svg.getBoundingClientRect()
      const w = rect.width || 600
      const h = rect.height || 200
      const canvas = document.createElement('canvas')
      canvas.width = w * scale
      canvas.height = h * scale
      const ctx = canvas.getContext('2d')
      ctx.fillStyle = '#fff'
      ctx.fillRect(0, 0, canvas.width, canvas.height)
      ctx.drawImage(img, 0, 0, canvas.width, canvas.height)
      URL.revokeObjectURL(url)
      canvas.toBlob((png) => {
        if (!png) return
        const a = document.createElement('a')
        a.href = URL.createObjectURL(png)
        a.download = 'motif.png'
        a.click()
        setTimeout(() => URL.revokeObjectURL(a.href), 1000)
      }, 'image/png')
    }
    img.src = url
  }

  return (
    <Stack spacing={1} alignItems="flex-start">
      <Box ref={containerRef} sx={{ width: '100%', overflowX: 'auto' }}>
        <DNALogo ppm={ppm} />
      </Box>
      <Button size="small" startIcon={<DownloadIcon />} onClick={downloadPNG}>
        Download PNG
      </Button>
    </Stack>
  )
}
