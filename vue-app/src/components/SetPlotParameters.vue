<!-- Based upon: https://vuematerial.io/components/form -->
<!-- sliders would be better for some inputs here, but the slider component is not yet implemented in this library -->
<!-- https://github.com/vuematerial/vue-material/blob/dev/ROADMAP.md -->

<template> 
  <div> 
    <md-progress-bar v-if="!initialPlotParametersSet || !modelParametersSet" md-mode="indeterminate" />
    <md-card 
      v-else
      style="margin: 10px auto;" 
    >
    
      <md-card-content>
        <div class="md-layout md-gutter">
          <div class="md-layout-item header-item" v-if="atLeastOneCanonicalTranscript">
            <div class="md-head">
              Canonical Transcript<span v-if="moreThanOneCanonicalTranscript">s</span>
            </div>
              <div v-for="canonicalTranscript in canonicalTranscripts" :key="canonicalTranscript.id"> 
                <md-divider></md-divider>
                <div class="md-subhead">
                  Name: {{ canonicalTranscript.display_name }}
                </div>
                <div class="md-subhead">
                  Ensemble ID: 
                  <a :href="getEnsembleUI(canonicalTranscript)" target="_blank">
                    {{ canonicalTranscript.id }}
                  </a>  
                </div>  
                <div class="md-subhead">
                  Strand: {{ canonicalTranscript.strand }}
                </div>
              </div>
          </div>

          <div class="md-layout-item header-item">
            <div class="md-head">
              Model 
            </div>
            <md-divider></md-divider>
            <div class="md-subhead">
              genomeBuild: {{ modelParameters.genomeBuild }} 
            </div>  
            <div class="md-subhead">
              kmerSize: {{ modelParameters.kmerSize }} 
            </div>  
            <div class="md-subhead">
              numberChromosomesMin: {{ modelParameters.numberChromosomesMin.toLocaleString() }} 
            </div>  
            <div class="md-subhead">
              windowSize: {{ modelParameters.windowSize }} 
            </div>  
          </div>

          <div class="md-layout-item header-item">
            <md-field>
              <label for="region">Region</label>
              <md-input
                id="region" 
                v-on:keyup.enter="getAPIData" 
                v-model="plotParameters.region" 
                :disabled="fetchingTimeSeriesData"
              />
            </md-field>
          </div>    

          <div class="md-layout-item" style="max-width: 150px;">
            <md-field>
              <label for="window-stride">Window stride</label>
              <md-input 
                id="window-stride" 
                v-on:keyup.enter="getAPIData" 
                v-model="plotParameters.windowStride" 
                :disabled="fetchingTimeSeriesData"
              />
            </md-field>
          </div> 
        </div>
        
      </md-card-content>

      <md-card-expand>
        <md-card-actions md-alignment="right">
          <md-button 
            v-on:click="getAPIData" 
            class="md-icon-button md-primary"
            :disabled="fetchingTimeSeriesData"
          >
            <md-icon>refresh</md-icon>
          </md-button>      

          <md-card-expand-trigger>
            <md-button  class="md-icon-button">
              <md-icon>keyboard_arrow_down</md-icon>
            </md-button>
          </md-card-expand-trigger>
        </md-card-actions>

        <md-card-expand-content>
          <md-card-content>
            <!-- <span v-if="atLeastOneCanonicalTranscript">Exons from the canonical transcript are labeled with their rank in that transcript.</span>
            <br> 
            Click on the time series plots to compute the null distribution of N and K at the corresponding genomic position.           
            <br>  -->
            <div class="md-head">
              DNA sequence
            </div>
            <md-divider></md-divider>
            <div 
              v-if="!fetchingTimeSeriesData && sequenceDataSet" 
              style="margin-top: 10px; width: 1500px; overflow-wrap: break-word;"
            >
              {{ sequenceData.sequence }}
            </div>
          </md-card-content>
        </md-card-expand-content>
      </md-card-expand>

      <md-snackbar 
        md-position="center" 
        :md-active.sync="showSnackbar" 
      >
        <span>Please correct the following error(s):</span>
        <ul>
          <li v-for="error in errors" :key="error">{{ error }}</li>
        </ul>    
        <md-button 
          class="md-primary" 
          @click="showSnackbar=false"
        >
          <md-icon>close</md-icon>
        </md-button>
      </md-snackbar>

    </md-card>
  </div>
</template> 

<script> 
// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Statements/import
import * as api from '@/api'

import { mapState, mapGetters } from 'vuex'

export default {
  name: 'SetPlotParameters',
  data () {
    return {
      errors: [],
      showSnackbar: false,

      plotParameters: null,
      initialPlotParametersSet: false,
    }
  }, 
  methods: {
    getEnsembleUI (canonicalTranscript) {
      return `https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=${canonicalTranscript.id}`
    },
    hasWhiteSpace (s) {
      return /\s/g.test(s)
    },
    getChromosomeStartEndFromWhiteSpaceDelimitedString (region) { 
      return region.split(/[ ]+/)
    },
    getChromosomeStartEndFromColonHyphenDelimitedString (region) { 
      const [chromosome, start_end] = region.split(':')
      const [start, end] = start_end.replaceAll(',', '').split('-')
      return [chromosome, start, end]
    },
    getChromosomeStartEnd (region) {
      if (this.hasWhiteSpace(region)) { 
        return this.getChromosomeStartEndFromWhiteSpaceDelimitedString(region)
      } else { 
        return this.getChromosomeStartEndFromColonHyphenDelimitedString(region)
      }
    },
    convertToColonHyphenDelimitedString (region) { 
      const [chromosome, start, end] = this.getChromosomeStartEnd(region)
      return `${chromosome}:${parseInt(start).toLocaleString()}-${parseInt(end).toLocaleString()}`
    },
    intervalSizeTooSmall (region) { 
      const [, start, end] = this.getChromosomeStartEnd(region)
      const intervalSize = parseInt(end) - parseInt(start)
      return intervalSize < this.modelParameters.windowSize
    },
    async getEndChromosomeLength (region) {
      const [chromosome, , end] = this.getChromosomeStartEnd(region)
      const chromosomeLength = await api.getChromosomeLength(chromosome, this.modelParameters.genomeBuild)
      return [parseInt(end), parseInt(chromosomeLength)]
    },
    async validateParameters () {
      this.plotParameters.region = this.convertToColonHyphenDelimitedString(this.plotParameters.region)

      this.errors = []

      if (!this.plotParameters.region) {
        this.errors.push('Region required.')
      }
      if (!this.plotParameters.windowStride) {
        this.errors.push('Window stride required.')
      }   

      if (this.intervalSizeTooSmall(this.plotParameters.region)) {
        this.errors.push(`Region must be at least ${this.modelParameters.windowSize}bp`)
      }

      const [end, chromosomeLength] = await this.getEndChromosomeLength(this.plotParameters.region)
      if ( end >= chromosomeLength ) {
        this.errors.push(`End coordinate, ${end}, must be smaller than chromosome length, ${chromosomeLength}`)
      }

      if (this.errors.length > 0) {
        return false
      } else { 
        return true
      }
    },
    async getAPIData () {
      if ( await this.validateParameters() ) {
        this.$store.dispatch('getExpectedObservedCounts', this.plotParameters)
        this.$store.dispatch('getCanonicalData', {
          region: this.plotParameters.region,
          genomeBuild: this.modelParameters.genomeBuild
        })
        this.$store.dispatch('getTrustworthyNoncodingRegions', this.plotParameters)
        this.$store.dispatch('getSequenceData', {
          region: this.plotParameters.region
        })
      } else { 
        this.showSnackbar = true
      }
    },
  },
  computed: {
    ...mapState([
      'canonicalTranscripts',
      'canonicalExons',
      'canonicalDataSet',
      'modelParameters',
      'modelParametersSet',
      'sequenceData',
      'sequenceDataSet'
    ]),
    ...mapGetters([
      'fetchingTimeSeriesData'
    ]),
    atLeastOneCanonicalTranscript () {
      return this.canonicalDataSet && Object.keys(this.canonicalTranscripts).length > 0
    },
    moreThanOneCanonicalTranscript () {
      return this.canonicalDataSet && Object.keys(this.canonicalTranscripts).length > 1
    }
  },
  async created () { 
    this.plotParameters = await api.getInitialPlotParameters()
    this.initialPlotParametersSet = true

    this.$store.dispatch('getModelParameters')
  }  
}
</script>

<style scoped> 
.header-item { 
  min-width: 250px; 
  max-width: 300px;
  padding-bottom: 10px;
}

div.md-subhead { 
  font-size: 12px;
}

.md-field.md-has-value .md-input { 
  font-size: 14px;
}
</style>

