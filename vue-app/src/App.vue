<template>
  <div style="padding: 10px;">
    <set-plot-parameters/>
    <plot-distributions/>
    <plot-time-series/>
  </div>
</template>

<script>

import PlotTimeSeries from '@/components/PlotTimeSeries'
import PlotDistributions from '@/components/PlotDistributions'
import SetPlotParameters from '@/components/SetPlotParameters'

export default {
  name: 'App',
  components: {
    PlotTimeSeries,
    PlotDistributions,
    SetPlotParameters
  }
}
</script>
