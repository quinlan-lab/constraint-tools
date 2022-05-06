<template>
  <div style="padding: 10px;">
    <set-plot-parameters/>
    <plot-individual-distributions/>
    <plot-time-series/>
  </div>
</template>

<script>

import PlotTimeSeries from '@/components/PlotTimeSeries'
import PlotIndividualDistributions from '@/components/PlotIndividualDistributions'
import SetPlotParameters from '@/components/SetPlotParameters'

export default {
  name: 'App',
  components: {
    PlotTimeSeries,
    PlotIndividualDistributions,
    SetPlotParameters
  }
}
</script>
