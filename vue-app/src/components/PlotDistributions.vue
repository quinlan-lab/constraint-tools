<template>
  <div>
    <!-- 
      Do not use "v-else" in the following element. 
      Otherwise the watcher callback will error
      because this element may not have yet rendered 
      when plotly code is called. 
      -->
    <div class="plot-container md-elevation-3">
      <div ref="plot"> </div>
    </div> 
    <div v-if="fetchingDistributionData" class="progress-bar-container">
      <md-progress-bar md-mode="indeterminate" />
    </div>
  </div>
</template>

<script>
import Plotly from 'plotly.js-dist'
import { mapState, mapGetters } from 'vuex'

export default {
  name: 'PlotTimeSeries',
  data () {
    return {
      yaxisMin: 1e-5, 
      yaxisMax: 1 
    }
  },
  methods: {
  },
  computed: {
    ...mapState([
      'NObserved',
      'selectedGenomicPosition',
      'distributionN',
      'modelParameters'
    ]),
    ...mapGetters([
      'fetchingDistributionData'
    ]),
    traces () {
      return [
        {
          x: [this.NObserved, this.NObserved],
          y: [this.yaxisMin, this.yaxisMax],
          name: 'observed # SNVs',
          xaxis: 'x',
          yaxis: 'y',
          mode: 'lines',
          line: {
            color: 'black',
            width: 1
          }
        },
        {
          x: this.distributionN['n'],
          y: this.distributionN['p(n)'],
          xaxis: 'x',
          yaxis: 'y',
          showlegend: false,
        }, 
        {
          x: [],
          y: [],
          name: 'Probability, p(k)',
          xaxis: 'x2',
          yaxis: 'y',
        }, 
      ] 
    },
    layout () {
      return { 
        showlegend: true,
        height: 300,
        width: 1000,
        grid: {
          rows: 1, 
          columns: 2, 
          subplots: [['xy'], ['x2y']],
        },
        xaxis: {
          domain: [0.0, 0.4],
          title: `# SNVs (per window), n`,
          showline: true,
          showgrid: false,
          zeroline: false,
          autotick: true,
          showticklabels: true,
        },
        xaxis2: {
          domain: [0.6, 1.0],
          title: ``,
          showline: true,
          showgrid: false,
          zeroline: false,
          autotick: true,
          showticklabels: true,
        },
        yaxis: { 
          title: `Null distribution, p(n), <br>at position ${this.selectedGenomicPosition.toLocaleString()}`,
          showgrid: false,
          showline: true,
          zeroline: false,
          autotick: true,
          showticklabels: true,
          type: 'log',
          // https://plotly.com/javascript/reference/layout/yaxis/#layout-yaxis-range
          range: [this.yaxisMin, this.yaxisMax].map(Math.log10),
          // showexponent: 'all',
          exponentformat: 'e',
        },
        responsive: true,
        font: {
          family: 'Roboto, sans-serif',
          size: 10
        },
        hovermode: 'closest',
        hoverlabel: {
          bgcolor: '#fafafa',
        },
        margin: { 
          t: 40 
        },
      }
    }
  },
  watch: {
    fetchingDistributionData: {
      handler: function (newValue, oldValue) {
        if ( newValue === true && oldValue === false ) { 
          console.log('fetching data for distribution plots from one or more APIs')
        }
        if ( newValue === false && oldValue === true ) {
          console.log('data for distribution plots fetched from all APIs')
          Plotly.react(this.$refs.plot, this.traces, this.layout)
        }
      },
      deep: true
    },
  },
  beforeUnmount () {
    Plotly.purge(this.$refs.plot)
  },
}
</script>

<style scoped>
  .plot-container { 
    margin: 10px auto; 
    background-color: white;
    position: relative;
  }

  .progress-bar-container {
    height: 10px;
  }
</style>
