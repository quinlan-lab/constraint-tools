<template>
  <div>
    <div v-if="fetchingAPIData" class="progress-bar-container">
      <md-progress-bar md-mode="indeterminate" />
    </div>
    <!-- 
      Do not use "v-else" in the following element. 
      Otherwise the watcher callback will error
      because this element may not have yet rendered 
      when plotly code is called. 
      -->
    <div class="plot-container md-elevation-3">
      <div ref="plot" style="padding-left: 60px"> </div>
      <div class="cpg-container">
        <div class="cpg-box" :style="{ 'background-color': cpgNegativeColor }"></div>
        <span class="cpg-explanation">CpG-</span>
        <br>
        <div class="cpg-box" :style="{ 'background-color': cpgPositiveColor }"></div>
        <span class="cpg-explanation">CpG+</span>
      </div> 
    </div> 
  </div>
</template>

<script>
import Plotly from 'plotly.js-dist'
import { mapState, mapGetters } from 'vuex'

export default {
  name: 'PlotCounts',
  data () {
    return {
      cpgNegativeColor: 'rgba(255, 0, 0, 0.2)',
      cpgPositiveColor: 'rgba(0, 0, 255, 0.2)'
    }
  },
  methods: {
    exonToRectangle (exon) {
      return {
        type: 'rect',
        xref: 'x', 
        yref: 'y2',
        x0: exon.start,
        y0: 0,
        x1: exon.end,
        y1: this.y2axisMax,
        fillcolor: '#d3d3d3',
        opacity: 0.6,
        line: {
            width: 0
        },
        layer: 'below',
      }
    },
    variantToEllipse (variant) {      
      return {
        type: 'circle',
        xref: 'x', 
        yref: 'y',
        x0: variant.position - this.semiMinorAxis,
        y0: -1 - Math.log10(variant.frequency),
        x1: variant.position + this.semiMinorAxis,
        y1: 1 + Math.log10(variant.frequency),
        fillcolor: variant.color, 
        line: {
            width: 0.5,
            color: 'rgba(0, 0, 0, 1)', 
        },
        // layer: 'below'
      }
    },
    ellipses (positions, frequencies, color) {
      const variants = positions.map((position, i) => (
        {
          'position': position, 
          'frequency': frequencies[i],
          'color': color
        }
      ))
      return variants.map(this.variantToEllipse)
    },
  },
  computed: {
    ...mapState([
      'expectedObservedCounts',
      'canonicalExons'
    ]),
    ...mapGetters([
      'fetchingAPIData'
    ]),
    semiMinorAxis () {
      const xaxisRange = this.expectedObservedCounts.end - this.expectedObservedCounts.start
      return 0.005 * xaxisRange
    },
    rectangles () {
      return this.canonicalExons.map(this.exonToRectangle) 
    },
    ellipsesCpGNegative () {
      const positions = this.expectedObservedCounts.snvCpGNegativePositions
      const frequencies = this.expectedObservedCounts.snvCpGNegativeFrequencies
      const color = this.cpgNegativeColor
      return this.ellipses(positions, frequencies, color)
    },
    ellipsesCpGPositive () {
      const positions = this.expectedObservedCounts.snvCpGPositivePositions
      const frequencies = this.expectedObservedCounts.snvCpGPositiveFrequencies
      const color = this.cpgPositiveColor
      return this.ellipses(positions, frequencies, color)
    },
    rectangleLabels () {
      return this.canonicalExons.map(exon => exon.rank) 
    },
    rectangleLabelXPositions () {
      return this.canonicalExons.map(exon => exon.start + 0.2*(exon.end - exon.start)) 
    },
    rectangleLabelYPositions () {
      return Array(this.canonicalExons.length).fill(0.95*this.y2axisMax)
    },
    y2axisMax () { 
      return Math.max(...this.expectedObservedCounts.NObserveds)
    },
    traces () {
      return [
        {
          x: this.expectedObservedCounts.windowPositions,
          y: this.expectedObservedCounts.NObserveds,
          name: 'N_observed',
          xaxis: 'x',
          yaxis: 'y2',
        }, 
        {
          x: this.rectangleLabelXPositions,
          y: this.rectangleLabelYPositions,
          text: this.rectangleLabels,
          mode: 'text',
          showlegend: false,
          xaxis: 'x',
          yaxis: 'y2'
        },      
        { 
          x: this.expectedObservedCounts.windowPositions,
          y: this.expectedObservedCounts.NBars,
          name: 'N_bar',
          xaxis: 'x',
          yaxis: 'y3',
        },
      ] 
    },
    layout () {
      return { 
        height: 750,
        grid: {
          rows: 3, 
          columns: 1, 
          subplots: [['xy'], ['xy2'], ['xy3']],
          roworder: 'top to bottom'
        },
        // https://codepen.io/plotly/pen/KpLVzv ...
        xaxis: {
          title: `Position (bps) along ${this.expectedObservedCounts.chromosome}`,
          showline: true,
          showgrid: false,
          zeroline: false,
          autotick: true,
          showticklabels: true,
          tickformat: ",.0f",
        },
        yaxis: { 
          domain: [0.8, 1.0],
          title: '1 + log(# ALT chromosomes) <br>(semi-major axis)',
          showgrid: false,
          showline: true,
          zeroline: false,
          autotick: true,
          showticklabels: true,
        },
        yaxis2: { 
          domain: [0.4, 0.75],
          title: 'observed #SNVs',
          showgrid: false,
          showline: true,
          zeroline: true,
          autotick: true,
          showticklabels: true,
          // range: [0, this.y2axisMax]
        },
        yaxis3: { 
          domain: [0, 0.35],
          title: 'z-scores',
          showgrid: false,
          showline: true,
          zeroline: false,
          autotick: true,
          showticklabels: true,
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
        shapes: [...this.rectangles, ...this.ellipsesCpGNegative, ...this.ellipsesCpGPositive]
      }
    }
  },
  watch: {
    fetchingAPIData: {
      handler: function (newValue, oldValue) {
        if ( newValue === true && oldValue === false ) { 
          console.log('fetching data from one or more APIs')
        }
        if ( newValue === false && oldValue === true ) {
          console.log('data fetched from all APIs')
          Plotly.react(this.$refs.plot, this.traces, this.layout)
        }
      },
      deep: true
    },
  },
  beforeUnmount () {
    Plotly.purge(this.$refs.plot)
  }  
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

  .cpg-container {
    position: absolute;
    top: 60px;
    left: 10px;
  }

  .cpg-box {
    width: 10px;
    height: 20px;
    border-radius: 50%;
    border: 1px solid black;
    margin: 5px;
    display: inline-block;
    vertical-align: middle;
  }

  .cpg-explanation {
    margin-left: 5px;
    vertical-align: middle;
  }

</style>
