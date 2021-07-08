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
    },
    traces () {
      return [
        {
          x: this.mutationCounts.windowPositions,
          y: this.mutationCounts.windowObservedMutationCounts,
          name: 'observed'
        }, 
        { 
          x: this.mutationCounts.windowPositions,
          y: this.mutationCounts.windowExpectedMutationCounts,
          name: 'expected'
        },
        {
          type: 'bar', // https://github.com/plotly/documentation/issues/1270#issuecomment-468645317
          x: this.mutationCounts.lollipopsCpGNegativePositions,
          y: this.mutationCounts.lollipopsCpGNegativeHeights,
          width: Array(this.mutationCounts.lollipopsCpGNegativePositions.length).fill(1),
          name: 'CpG-'
        },
        {
          type: 'bar', // https://github.com/plotly/documentation/issues/1270#issuecomment-468645317
          x: this.mutationCounts.lollipopsCpGPositivePositions,
          y: this.mutationCounts.lollipopsCpGPositiveHeights,
          width: Array(this.mutationCounts.lollipopsCpGPositivePositions.length).fill(1),
          name: 'CpG+'
        },
      ] 
    },
    layout () {
      return { 
        responsive: true,
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
      }
    }
  },
  watch: {
    mutationCounts: {
      handler: function () {
        Plotly.react(this.$refs.chart, this.traces, this.layout)
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
  margin: 10px auto; 
  width: 1000px; 
  height: 500px;
  background-color: white;
}

.progress-bar-container {
  height: 10px;
}
</style>
