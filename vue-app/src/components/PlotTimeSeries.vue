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
        <span class="annotation"># chroms > {{ modelParameters.numberChromosomesMin.toLocaleString() }}</span>
        <br>
        <div class="cpg-box" :style="{ 'background-color': cpgNegativeColor }"></div>
        <span class="explanation">CpG-</span>
        <br>
        <div class="cpg-box" :style="{ 'background-color': cpgPositiveColor }"></div>
        <span class="explanation">CpG+</span>
        <div style="padding-top: 200px">
          <div class="box" :style="{ 'background-color': exonColor }"></div>
          <span class="explanation">exon</span>
          <br>
          <div class="box" :style="{ 'background-color': trustworthyNoncodingRegionColor }"></div>
          <span class="explanation">trustworthy noncoding region</span>
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
      y2axisMin: 0,
      y3axisMin: -10,
      y3axisMax: 10,
      plotlyEventListenersAdded: false,
      viewBreakpoint: 150
    }
  },
  methods: {
    getLayoutShowSequence () {
      return { ...this.layoutCore, xaxis: this.xaxisShowSequence }
    },
    getLayoutShowCoordinates () {
      // converting this to a method results in a bug:
      // the "xaxis" property sporadically changes to the value of "this.xaxisShowSequence" 
      return { ...this.layoutCore, xaxis: this.xaxisShowCoordinates }
    },
    relayout (start, end) { 
      const length = end - start
      // https://github.com/plotly/plotly.js/issues/1010#issuecomment-252952758
      Plotly.relayout(
        this.$refs.plot, 
        length < this.viewBreakpoint ? this.getLayoutShowSequence() : this.getLayoutShowCoordinates()
      )
    },
    hasProperty (object, property) {
      return Object.prototype.hasOwnProperty.call(object, property)
    },
    handleRelayoutEvent (eventData) {
      // event-data definition: https://plotly.com/javascript/plotlyjs-events/#update-data
      console.log('relayout-event data:')
      console.log(eventData)

      if (this.hasProperty(eventData, 'xaxis.autorange')) {
        this.relayout(this.xaxisMin, this.xaxisMax)
        return
      }

      if (
        this.hasProperty(eventData, 'xaxis.range[0]') &&
        this.hasProperty(eventData, 'xaxis.range[1]')
      ) {
        const start = parseInt(eventData['xaxis.range[0]'])
        const end = parseInt(eventData['xaxis.range[1]'])
        this.relayout(start, end)
        return
      }
    },
    getDistributions (eventData) {
      const [point] = eventData.points
      this.$store.commit('setSelectedGenomicPosition', point.x)
      const index = point.pointNumber

      this.$store.dispatch('getDistributions', {
        'region': this.expectedObservedCounts.windowRegions[index],
        'M': this.expectedObservedCounts.Ms[index]
      })

      this.$store.commit('setObservedValues', { 
        'N': this.expectedObservedCounts.NObserveds[index],
        'K': this.expectedObservedCounts.KObserveds[index]
      })
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
    computeMeanStdSNVCount () {
      const n = this.modelParameters.windowSize // number of Bernoulli trials, each being an attempt to substitute a nucleotide 
      const p = 0.3 // probability of success: typical substitution probability
      const mean = n*p // (binomial) mean SNV count
      const std = Math.sqrt(n*p*(1-p)) // (binomial) std of SNV count
      return [mean, std]
    },
    y2axisMax () {
      const [mean, std] = this.computeMeanStdSNVCount
      console.log(`predicted mean SNV count: ${mean}`)
      console.log(`predicted std SNV count: ${std}`)
      return mean + 5*std
    },
    nucleotides () {
      return [...this.sequenceData.sequence]
    },
    nucleotidePositions () {
      const xs = []
      for (let x = this.sequenceData.start; x < this.sequenceData.end; x++) {
        xs.push(x)
      }      
      return xs
    },
    fetchingAnyData () {
      return this.fetchingTimeSeriesData || this.fetchingDistributions
    },
    ...mapState([
      'expectedObservedCounts',
      'canonicalExons',
      'modelParameters',
      'trustworthyNoncodingRegions',
      'sequenceData',
      'selectedGenomicPosition',
      'fetchingDistributions',
      'exonColor',
      'trustworthyNoncodingRegionColor'
    ]),
    ...mapGetters([
      'fetchingTimeSeriesData'
    ]),
    semiMinorAxis () {
      const xaxisRange = this.expectedObservedCounts.end - this.expectedObservedCounts.start
      return 0.005 * xaxisRange
    },
    exonRectangles () {
      return this.canonicalExons.map(exon => this.createRectangle(exon, this.exonColor)) 
    },
    trustworthyNoncodingRegionRectangles () {
      return this.trustworthyNoncodingRegions.map(
        trustworthyNoncodingRegion => this.createRectangle(
          trustworthyNoncodingRegion, 
          this.trustworthyNoncodingRegionColor
      ))
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
          name: 'position at which null distributions of N and K are evaluated',
          xaxis: 'x',
          yaxis: 'y3',
          mode: 'lines',
          line: {
            color: 'black',
            width: 2
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
            width: 2
          }
        },
        {
          x: this.expectedObservedCounts.windowPositions,
          y: this.expectedObservedCounts.NObserveds,
          name: 'N_observed (no constraint on chromosome number)',
          xaxis: 'x',
          yaxis: 'y2',
          line: {
            width: 1
          }
        }, 
        {
          x: this.expectedObservedCounts.windowPositions,
          y: this.expectedObservedCounts.Ms,
          name: `M (# chroms > ${this.modelParameters.numberChromosomesMin.toLocaleString()})`,
          xaxis: 'x',
          yaxis: 'y2',
          line: {
            width: 1
          }
        }, 
        {
          x: this.expectedObservedCounts.windowPositions,
          y: this.expectedObservedCounts.KObserveds,
          name: `singletons (K_observed; # chroms > ${this.modelParameters.numberChromosomesMin.toLocaleString()})`,
          xaxis: 'x',
          yaxis: 'y2',
          line: {
            width: 1
          }
        }, 
        { 
          x: this.expectedObservedCounts.windowPositions,
          y: this.expectedObservedCounts.NBars,
          name: 'N_bar (no constraint on chromosome number)',
          xaxis: 'x',
          yaxis: 'y3',
          line: {
            width: 1
          }
        },
        { 
          x: this.expectedObservedCounts.windowPositions,
          y: this.expectedObservedCounts.KBars,
          name: `K_bar (# chroms > ${this.modelParameters.numberChromosomesMin.toLocaleString()})`,
          xaxis: 'x',
          yaxis: 'y3',
          line: {
            width: 1
          }
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
    xaxisCore () {
      return {
        title: `Position (bps) along ${this.expectedObservedCounts.chromosome}`,
        showline: true,
        showgrid: false,
        zeroline: false,
        showticklabels: true,
        range: [this.xaxisMin, this.xaxisMax],
      }
    },
    xaxisShowSequence () {
      return {
        ...this.xaxisCore,  
        tickmode: 'array',
        ticktext: this.nucleotides,
        tickvals: this.nucleotidePositions,
        tickangle: 0,
      }
    },
    xaxisShowCoordinates () {
      return {
        ...this.xaxisCore,  
        tickmode: 'auto',
        tickformat: ",.0f",
        autotick: true,
      }
    },
    layoutCore () {
      return { 
        showlegend: true,
        legend: {
          x: 1.05,
          y: 0.4,
          // xanchor: 'right',
        },
        height: 500,
        width: 1400,
        grid: {
          rows: 3, 
          columns: 1, 
          subplots: [['xy'], ['xy2'], ['xy3']],
          roworder: 'top to bottom'
        },
        // https://codepen.io/plotly/pen/KpLVzv ...
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
          title: 'obs. # SNVs <br>(per window)',
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
          zeroline: true,
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
        shapes: [
          ...this.exonRectangles, 
          ...this.trustworthyNoncodingRegionRectangles, 
          ...this.cpgNegativeEllipses, 
          ...this.cpgPositiveEllipses
        ]
      }
    },
  },
  watch: {
    fetchingAnyData: {
      handler: function (newValue, oldValue) {
        if ( newValue === true && oldValue === false ) { 
          console.log('fetching data from one or more APIs')
        }
        if ( newValue === false && oldValue === true ) {
          console.log('data fetched from all APIs')

          const xaxisLength = this.xaxisMax - this.xaxisMin
          Plotly.react(
            this.$refs.plot, 
            this.traces, 
            xaxisLength < this.viewBreakpoint ? this.getLayoutShowSequence() : this.getLayoutShowCoordinates()
          )

          if ( !this.plotlyEventListenersAdded ) { 
            this.$refs.plot.on('plotly_click', this.getDistributions)
            this.$refs.plot.on('plotly_relayout', this.handleRelayoutEvent)
            console.log('plotly event listeners added')
            this.plotlyEventListenersAdded = true
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
    top: 30px;
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
