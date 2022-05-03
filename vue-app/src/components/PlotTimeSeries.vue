<template>
  <div>
    <div v-if="fetchingTimeSeriesData" class="progress-bar-container">
      <md-progress-bar md-mode="indeterminate" />
    </div>
    <!-- 
      Do not use "v-else" in the following element. 
      Otherwise the watcher callback will error
      because this element may not have yet rendered 
      when plotly code is called. 
      -->
    <div class="plot-container md-elevation-3">
      <div ref="plot" style="padding-left: 120px"> </div>
      <div class="legend-container" v-if="expectedObservedCounts">
        <span class="annotation"># chroms > {{ modelParameters.numberChromosomesMin }}</span>
        <br>
        <div class="cpg-box" :style="{ 'background-color': cpgNegativeColor }"></div>
        <span class="explanation">CpG-</span>
        <br>
        <div class="cpg-box" :style="{ 'background-color': cpgPositiveColor }"></div>
        <span class="explanation">CpG+</span>
        <div style="padding-top: 250px">
          <div class="box" :style="{ 'background-color': exonColor }"></div>
          <span class="explanation">exon</span>
          <br>
          <div class="box" :style="{ 'background-color': neutralRegionColor }"></div>
          <span class="explanation">"neutral region"</span>
        </div>
      </div> 
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
      cpgNegativeColor: 'rgba(255, 0, 0, 0.2)',
      cpgPositiveColor: 'rgba(0, 0, 255, 0.2)',
      exonColor: '#d3d3d3',
      neutralRegionColor: 'rgba(0, 255, 0, 0.3)',
      y2axisMin: 0,
      y2axisMax: 50,
      y3axisMin: -10,
      y3axisMax: 10,
      plotlyEventListenerAdded: false
    }
  },
  methods: {
    getDistributionN (eventData) {
      const [point] = eventData.points
      this.$store.commit('setSelectedGenomicPosition', point.x)
      const index = point.pointNumber

      const window = { 
        position: this.expectedObservedCounts.windowPositions[index], 
        region: this.expectedObservedCounts.windowRegions[index]
      }
      console.log('window:')
      console.log(window)
      this.$store.dispatch('getDistributionN', window)

      const NObserved = this.expectedObservedCounts.NObserveds[index]
      console.log('NObserved:')
      console.log(NObserved)
      this.$store.commit('setNObserved', NObserved)
    },
    createRectangle (region, color) {
      return {
        type: 'rect',
        xref: 'x', 
        yref: 'y3',
        x0: region.start,
        y0: this.y3axisMin,
        x1: region.end,
        y1: this.y3axisMax,
        fillcolor: color,
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
    fetchingAnyData () {
      return this.fetchingTimeSeriesData || this.fetchingDistributionData
    },
    ...mapState([
      'expectedObservedCounts',
      'canonicalExons',
      'modelParameters',
      'neutralRegions',
      'selectedGenomicPosition'
    ]),
    ...mapGetters([
      'fetchingTimeSeriesData',
      'fetchingDistributionData'
    ]),
    semiMinorAxis () {
      const xaxisRange = this.expectedObservedCounts.end - this.expectedObservedCounts.start
      return 0.005 * xaxisRange
    },
    exonRectangles () {
      return this.canonicalExons.map(exon => this.createRectangle(exon, this.exonColor)) 
    },
    neutralRegionRectangles () {
      return this.neutralRegions.map(neutralRegion => this.createRectangle(neutralRegion, this.neutralRegionColor))
    },
    cpgNegativeEllipses () {
      const positions = this.expectedObservedCounts.snvCpGNegativePositions
      const frequencies = this.expectedObservedCounts.snvCpGNegativeFrequencies
      const color = this.cpgNegativeColor
      return this.ellipses(positions, frequencies, color)
    },
    cpgPositiveEllipses () {
      const positions = this.expectedObservedCounts.snvCpGPositivePositions
      const frequencies = this.expectedObservedCounts.snvCpGPositiveFrequencies
      const color = this.cpgPositiveColor
      return this.ellipses(positions, frequencies, color)
    },
    exonRectangleLabels () {
      return this.canonicalExons.map(exon => exon.rank) 
    },
    exonRectangleLabelXPositions () {
      return this.canonicalExons.map(exon => exon.start + 0.2*(exon.end - exon.start)) 
    },
    exonRectangleLabelYPositions () {
      return Array(this.canonicalExons.length).fill(0.9*this.y3axisMax)
    },
    xaxisMin () { 
      return Math.min(...this.expectedObservedCounts.windowPositions)
    },
    xaxisMax () { 
      return Math.max(...this.expectedObservedCounts.windowPositions)
    },
    traces () {
      return [
        {
          x: [this.selectedGenomicPosition, this.selectedGenomicPosition],
          y: [this.y3axisMin, this.y3axisMax],
          name: 'position at which null distributions <br>of N and K is evaluated',
          xaxis: 'x',
          yaxis: 'y3',
          mode: 'lines',
          line: {
            color: 'black',
            width: 1
          }
        },
        {
          x: [this.selectedGenomicPosition, this.selectedGenomicPosition],
          y: [this.y2axisMin, this.y2axisMax],
          showlegend: false,
          xaxis: 'x',
          yaxis: 'y2',
          mode: 'lines',
          line: {
            color: 'black',
            width: 1
          }
        },
        {
          x: this.expectedObservedCounts.windowPositions,
          y: this.expectedObservedCounts.NObserveds,
          name: 'N_observed (no constraint <br> on chromosome number)',
          xaxis: 'x',
          yaxis: 'y2',
        }, 
        { 
          x: this.expectedObservedCounts.windowPositions,
          y: this.expectedObservedCounts.NBars,
          name: 'N_bar (no constraint <br> on chromosome number)',
          xaxis: 'x',
          yaxis: 'y3',
        },
        {
          x: this.exonRectangleLabelXPositions,
          y: this.exonRectangleLabelYPositions,
          text: this.exonRectangleLabels,
          mode: 'text',
          showlegend: false,
          xaxis: 'x',
          yaxis: 'y3'
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
        height: 600,
        width: 1000,
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
          range: [this.xaxisMin, this.xaxisMax]
        },
        yaxis: { 
          domain: [0.8, 1.0],
          title: '1 + log(# ALT chroms) <br>(semi-major axis)',
          showgrid: false,
          showline: true,
          zeroline: false,
          autotick: true,
          showticklabels: true,
        },
        yaxis2: { 
          domain: [0.4, 0.75],
          title: 'observed # SNVs <br>(per window)',
          showgrid: false,
          showline: true,
          zeroline: true,
          autotick: true,
          showticklabels: true,
          range: [this.y2axisMin, this.y2axisMax]
        },
        yaxis3: { 
          domain: [0, 0.35],
          title: 'z-scores',
          showgrid: false,
          showline: true,
          zeroline: false,
          autotick: true,
          showticklabels: true,
          range: [this.y3axisMin, this.y3axisMax]
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
        shapes: [...this.exonRectangles, ...this.neutralRegionRectangles, ...this.cpgNegativeEllipses, ...this.cpgPositiveEllipses]
      }
    }
  },
  watch: {
    fetchingAnyData: {
      handler: function (newValue, oldValue) {
        if ( newValue === true && oldValue === false ) { 
          console.log('fetching data from one or more APIs')
        }
        if ( newValue === false && oldValue === true ) {
          console.log('data fetched from all APIs')
          Plotly.react(this.$refs.plot, this.traces, this.layout)
          if ( !this.plotlyEventListenerAdded ) { 
            console.log('plotly event listener added')
            this.$refs.plot.on('plotly_click', this.getDistributionN)
            this.plotlyEventListenerAdded = true
          }
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

  .legend-container {
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

  .box {
    width: 10px;
    height: 20px;
    border-radius: 0%;
    border: 0px solid black;
    margin: 5px;
    display: inline-block;
    vertical-align: middle;
  }

  .explanation {
    margin-left: 5px;
    vertical-align: middle;
  }

  .annotation {
    color: grey;
  }
</style>
