<template>
  <div>
    <set-plot-parameters/>
    <plot-mutation-counts/>
  </div>
</template>

<script>

import PlotMutationCounts from '@/components/PlotMutationCounts'
import SetPlotParameters from '@/components/SetPlotParameters'

export default {
  name: 'App',
  components: {
    PlotMutationCounts,
    SetPlotParameters
  }
}
</script>
