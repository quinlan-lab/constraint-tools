<template>
  <div class="plot-container md-elevation-3">
    <div class="progress-bar-container">
      <md-progress-bar v-if="fetchingData" md-mode="indeterminate" />
    </div>
    <div ref="chart"></div>
  </div> 
</template>

<script>
import Plotly from 'plotly.js-dist'
import { mapState } from 'vuex'

export default {
  name: 'PlotMutationCounts',
  computed: {
    ...mapState([
      'mutationCounts',
      'fetchingData'
    ]),
    yaxisMax () { 
      return Math.max(
        ...this.mutationCounts.windowExpectedMutationCounts,
        ...this.mutationCounts.windowObservedMutationCounts
      )
    }
  },
  watch: {
    mutationCounts: {
      handler: function () {
        Plotly.react(
          this.$refs.chart,
          [
            {
              x: this.mutationCounts.windowPositions,
              y: this.mutationCounts.windowObservedMutationCounts,
              name: 'observed'
            }, 
            { 
              x: this.mutationCounts.windowPositions,
              y: this.mutationCounts.windowExpectedMutationCounts,
              name: 'expected'
            }
          ], 
          { 
            font: {
              family: 'Roboto, sans-serif',
              size: 16
            },
            // https://codepen.io/plotly/pen/KpLVzv ...
            xaxis: {
              title: `Position (bps) along chromosome ${this.mutationCounts.chromosome}`,
              showline: true,
              showgrid: false,
              zeroline: false,
              autotick: true,
              showticklabels: true,
            },
            yaxis: {
              title: 'Number of mutations',
              showline: true,
              showgrid: false,
              zeroline: false,
              autotick: true,
              showticklabels: true,
              range: [0, this.yaxisMax]
            },
            hovermode: 'closest',
            hoverlabel: {
              bgcolor: '#fafafa'
            },
            margin: { 
              t: 40 
            } 
          },
        )
      },
      deep: true
    }
  },
  beforeUnmount () {
    Plotly.purge(this.$refs.chart)
  }  
}
</script>

<style scoped>
.plot-container { 
  margin: auto; 
  width: 800px; 
  height: 500px;
  background-color: white;
}

.progress-bar-container {
  height: 10px;
}
</style>
