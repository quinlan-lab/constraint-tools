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
    <div v-if="fetchingDistributions" class="progress-bar-container">
      <md-progress-bar md-mode="indeterminate" />
    </div>
  </div>
</template>

<script>
import Plotly from 'plotly.js-dist'
import { mapState } from 'vuex'

export default {
  name: 'PlotIndividualDistributions',
  data () {
    return {
    }
  },
  methods: {
  },
  computed: {
    ...mapState([
      'observedValues',
      'selectedGenomicPosition',
      'distributions',
      'modelParameters',
      'fetchingDistributions'
    ]),
    distributionN () {
      return this.distributions['N']
    },
    distributionK () {
      return this.distributions['K']
    },
    yaxisMax () {
      return Math.max(...this.distributionN['p(n)'])
    },
    yaxis2Max () {
      return Math.max(...this.distributionK['p(k)'])
    },
    traces () {
      return [
        {
          x: [this.observedValues.N, this.observedValues.N],
          y: [0, this.yaxisMax],
          name: 'observed values',
          xaxis: 'x',
          yaxis: 'y',
          mode: 'lines',
          line: {
            color: 'black',
            width: 2
          },
        },
        {
          x: this.distributionN['n'],
          y: this.distributionN['p(n)'],
          mode: 'lines',
          xaxis: 'x',
          yaxis: 'y',
          showlegend: false,
          line: {
            color: 'rgba(0,0,0,0.3)',
            width: 1
          }
        }, 
        {
          x: [this.observedValues.K, this.observedValues.K],
          y: [0, this.yaxis2Max],
          name: '',
          showlegend: false,
          xaxis: 'x2',
          yaxis: 'y2',
          mode: 'lines',
          line: {
            color: 'black',
            width: 2
          }
        },
        {
          x: this.distributionK['k'],
          y: this.distributionK['p(k)'],
          mode: 'lines',
          xaxis: 'x2',
          yaxis: 'y2',
          showlegend: false,
          line: {
            color: 'rgba(0,0,0,0.3)',
            width: 1
          }
        }, 
      ] 
    },
    layout () {
      return { 
        annotations: [
          {
            xref: 'x',
            yref: 'y',
            x: this.observedValues.N,
            xanchor: 'middle',
            y: this.yaxisMax,
            yanchor: 'bottom',
            text: this.observedValues.N,
            showarrow: false
          },
          {
            xref: 'x2',
            yref: 'y2',
            x: this.observedValues.K,
            xanchor: 'middle',
            y: this.yaxis2Max,
            yanchor: 'bottom',
            text: this.observedValues.K,
            showarrow: false
          },
          {
            xref: 'x2',
            yref: 'y2',
            x: 0,
            xanchor: 'left',
            y: 1.2*this.yaxis2Max, 
            yanchor: 'bottom',
            text: `# windows = ${this.distributionK['windowCount'].toLocaleString()}`,
            showarrow: false
          }
        ],
        showlegend: true,
        height: 250,
        width: 1000,
        grid: {
          rows: 1, 
          columns: 2, 
          pattern: 'independent',
        },
        xaxis: {
          domain: [0.0, 0.4],
          title: `# SNVs (per window), n`,
          showline: true,
          showgrid: false,
          zeroline: false,
          autotick: true,
          showticklabels: true,
          rangemode : 'tozero'
        },
        xaxis2: {
          domain: [0.6, 1.0],
          title: `# singletons (per window), k`,
          showline: true,
          showgrid: false,
          zeroline: false,
          autotick: true,
          showticklabels: true,
          rangemode : 'tozero'
        },
        yaxis: { 
          title: `Null distribution, p(n), <br>at position ${this.selectedGenomicPosition.toLocaleString()}`,
          showgrid: false,
          showline: true,
          zeroline: false,
          autotick: true,
          showticklabels: true,
          rangemode : 'tozero'
          // type: 'log',
          // // https://plotly.com/javascript/reference/layout/yaxis/#layout-yaxis-range
          // range: [this.yaxisMin, this.yaxisMax].map(Math.log10),
          // // showexponent: 'all',
          // exponentformat: 'e',
        },
        yaxis2: { 
          title: `Null distribution, p(k), <br>at position ${this.selectedGenomicPosition.toLocaleString()}`,
          showgrid: false,
          showline: true,
          zeroline: false,
          autotick: true,
          showticklabels: true,
          rangemode : 'tozero'
          // type: 'log',
          // // https://plotly.com/javascript/reference/layout/yaxis/#layout-yaxis-range
          // range: [this.yaxisMin, this.yaxisMax].map(Math.log10),
          // // showexponent: 'all',
          // exponentformat: 'e',
        },
        responsive: true,
        font: {
          family: 'Roboto, sans-serif',
          size: 12
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
    fetchingDistributions: {
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
