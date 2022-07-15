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
      'exonColor',
      'trustworthyNoncodingRegionColor'
    ]),
    ...mapGetters([
      'fetchingTimeSeriesData'
    ]),
    traces () {
      return [
        { 
          x: this.expectedObservedCounts.NBarsTrustworthyNoncodingRegions,
          y: this.expectedObservedCounts.KBarsTrustworthyNoncodingRegions,
          name: 'trustworthy noncoding region',
          xaxis: 'x',
          yaxis: 'y',
          mode: 'markers',
          type: 'scatter',
          marker: { 
            size: 5,
            color: this.trustworthyNoncodingRegionColor
          },
          showlegend: true,
        },
        { 
          x: this.expectedObservedCounts.NBarsExons,
          y: this.expectedObservedCounts.KBarsExons,
          name: 'exon',
          xaxis: 'x',
          yaxis: 'y',
          mode: 'markers',
          type: 'scatter',
          marker: { 
            size: 5,
            color: this.exonColor
          },
          showlegend: true,
        },
      ] 
    },
    layout () {
      return { 
        showlegend: true,
        legend: {
          x: 1.05,
          y: 0.4,
          // xanchor: 'right',
        },
        height: 300,
        width: 500,
        xaxis: {
          title: `Z-score for <br># SNVs (per window), <br>N_bar`,
          showline: false,
          showgrid: false,
          zeroline: true,
          autotick: true,
          showticklabels: true,
          range: [-10, 10]
        },
        yaxis: { 
          title: 'Z-score for <br># singletons (per window), <br>K_bar',
          showgrid: false,
          showline: false,
          zeroline: true,
          autotick: true,
          showticklabels: true,
          range: [-5, 5],
          // scaleanchor: "x",
          // scaleratio: 1,

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
  }

</style>
