<template>
  <div style="padding: 10px;">
    <set-plot-parameters />
    <div class="md-layout md-gutter">
      <plot-individual-distributions class="md-layout-item"/>
      <plot-joint-distribution class="md-layout-item"/>
    </div>
    <plot-time-series/>
  </div>
</template>

<script>

import PlotTimeSeries from '@/components/PlotTimeSeries'
import PlotIndividualDistributions from '@/components/PlotIndividualDistributions'
import PlotJointDistribution from '@/components/PlotJointDistribution'
import SetPlotParameters from '@/components/SetPlotParameters'

export default {
  name: 'App',
  components: {
    PlotTimeSeries,
    PlotIndividualDistributions,
    PlotJointDistribution,
    SetPlotParameters
  }
}
</script>

<style scoped>
</style>
