<template>
  <div>
    <div class="plot-container md-elevation-3">
      <div ref="plot"> </div>
    </div> 
  </div>
</template>

<script>
import Plotly from 'plotly.js-dist'
import { mapState, mapGetters } from 'vuex'

export default {
  name: 'PlotJointDistribution',
  data () {
    return {
    }
  },
  methods: {
  },
  computed: {
    ...mapState([
      'expectedObservedCounts',
    ]),
    ...mapGetters([
      'fetchingTimeSeriesData'
    ]),
    traces () {
      return [
        { 
          x: this.expectedObservedCounts.NBars,
          y: this.expectedObservedCounts.KBars,
          name: '',
          xaxis: 'x',
          yaxis: 'y',
          mode: 'markers',
          type: 'scatter',
          marker: { 
            size: 5,
            color: 'black'
          }
        },
      ] 
    },
    layout () {
      return { 
        showlegend: false,
        legend: {
          x: 1.05,
          y: 0.4,
          // xanchor: 'right',
        },
        height: 300,
        width: 300,
        xaxis: {
          title: `Z-score for <br># SNVs (per window), <br>N_bar`,
          showline: true,
          showgrid: false,
          zeroline: false,
          autotick: true,
          showticklabels: true,
          // range: [-5, 5]
        },
        yaxis: { 
          title: 'Z-score for <br># singletons (per window), <br>K_bar',
          showgrid: false,
          showline: true,
          zeroline: false,
          autotick: true,
          showticklabels: true,
          // range: [-2, 2]
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
          t: 40,
        },
      }
    }
  },
  watch: {
    fetchingTimeSeriesData: {
      handler: function (newValue, oldValue) {
        if ( newValue === true && oldValue === false ) { 
          console.log('fetching time-series data from one or more APIs')
        }
        if ( newValue === false && oldValue === true ) {
          console.log('time-series data fetched from all APIs')
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
