<template>
  <div class="plot-container md-elevation-3">
    <div v-if="fetchingAPIData" class="progress-bar-container">
      <md-progress-bar md-mode="indeterminate" />
    </div>
    <!-- 
      Do not use "v-else" in the following element. 
      Otherwise the watcher callback will error
      because this element may not have yet rendered 
      when plotly code is called. 
      -->
    <div ref="chart"> </div>
  </div> 
</template>

<script>
import Plotly from 'plotly.js-dist'
import { mapState, mapGetters } from 'vuex'

export default {
  name: 'PlotMutationCounts',
  methods: {
    exonToRectangle (exon) {
      return {
        type: 'rect',
        xref: 'x', 
        yref: 'y',
        x0: exon.start,
        y0: 0,
        x1: exon.end,
        y1: this.yaxisMax,
        fillcolor: '#d3d3d3',
        opacity: 0.6,
        line: {
            width: 0
        },
        layer: 'below'
      }
    } 
  },
  computed: {
    ...mapState([
      'mutationCounts',
      'canonicalExons'
    ]),
    ...mapGetters([
      'fetchingAPIData'
    ]),
    rectangles () {
      return this.canonicalExons.map(this.exonToRectangle) 
    },
    rectangleLabels () {
      return this.canonicalExons.map(exon => exon.rank) 
    },
    rectangleLabelXPositions () {
      return this.canonicalExons.map(exon => exon.start + 0.2*(exon.end - exon.start)) 
    },
    rectangleLabelYPositions () {
      return Array(this.canonicalExons.length).fill(0.95*this.yaxisMax)
    },
    yaxisMax () { 
      return Math.max(
        ...this.mutationCounts.windowExpectedMutationCounts,
        ...this.mutationCounts.windowObservedMutationCounts
      )
    },
    traces () {
      return [
        {
          x: this.rectangleLabelXPositions,
          y: this.rectangleLabelYPositions,
          text: this.rectangleLabels,
          mode: 'text',
          showlegend: false
        },      
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
        },
        shapes: this.rectangles
      }
    }
  },
  watch: {
    fetchingAPIData: {
      handler: function (newValue, oldValue) {
        if ( newValue === true && oldValue === false) { 
          console.log('fetching data from one or more APIs')
        }
        if ( newValue === false && oldValue === true ) {
          console.log('data fetched from all APIs')
          Plotly.react(this.$refs.chart, this.traces, this.layout)
        }
      },
      deep: true
    },
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
