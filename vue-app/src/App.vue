<template>
  <div style="padding: 10px;">
    <set-plot-parameters />
    <div class="md-layout">
      <plot-individual-distributions class="md-layout-item distribution"/>
      <plot-joint-distribution class="md-layout-item distribution"/>
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
  .distribution { 
    margin-left: 10px; 
    margin-right: 10px; 
  }
</style>
