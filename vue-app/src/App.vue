<template>
  <div style="padding: 10px;">
    <set-plot-parameters/>
    <plot-counts/>
  </div>
</template>

<script>

import PlotCounts from '@/components/PlotCounts'
import SetPlotParameters from '@/components/SetPlotParameters'

export default {
  name: 'App',
  components: {
    PlotCounts,
    SetPlotParameters
  }
}
</script>
