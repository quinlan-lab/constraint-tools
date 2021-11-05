<!-- Based upon: https://vuematerial.io/components/form -->
<!-- sliders would be better for some inputs here, but the slider component is not yet implemented in this library -->
<!-- https://github.com/vuematerial/vue-material/blob/dev/ROADMAP.md -->

<template> 
  <div> 
    <md-progress-bar v-if="!initialPlotParametersSet" md-mode="indeterminate" />
    <md-card 
      v-else
      style="margin: 10px auto;" 
    >
    
      <md-card-header v-if="canonicalTranscriptExists"> 
        <div class="md-title">
          Canonical Transcript: {{ canonicalTranscript.display_name }}
        </div>
        <div class="md-subhead">
          Ensemble ID: 
          <a :href="canonicalTranscriptEnsembleUI" target="_blank">
            {{ canonicalTranscript.id }}
          </a>  
        </div>  
        <div class="md-subhead">
          Strand: {{ canonicalTranscript.strand }}
        </div>
      </md-card-header>

      <md-card-content>
        <div class="md-layout md-gutter">
          <div class="md-layout-item" style="min-width: 250px; max-width: 300px;">
            <md-field>
              <label for="region">Region</label>
              <md-input id="region" v-model="plotParameters.region" :disabled="fetchingAPIData" />
            </md-field>
          </div>    

          <div class="md-layout-item" style="max-width: 150px;">
            <md-field>
              <label for="window-size">Window size</label>
              <md-input id="window-size" v-model="plotParameters.windowSize" :disabled="fetchingAPIData" />
            </md-field>
          </div>    

          <div class="md-layout-item" style="max-width: 150px;">
            <md-field>
              <label for="window-stride">Window stride</label>
              <md-input id="window-stride" v-model="plotParameters.windowStride" :disabled="fetchingAPIData" />
            </md-field>
          </div> 
        </div>
        
      </md-card-content>

      <md-card-expand>
        <md-card-actions md-alignment="right">
          <md-button 
            v-on:click="getAPIData" 
            class="md-icon-button md-primary"
            :disabled="fetchingAPIData"
          >
            <md-icon>refresh</md-icon>
          </md-button>      

          <md-card-expand-trigger v-if="canonicalTranscriptExists">
            <md-button  class="md-icon-button">
              <md-icon>keyboard_arrow_down</md-icon>
            </md-button>
          </md-card-expand-trigger>
        </md-card-actions>

        <md-card-expand-content>
          <md-card-content>
            Exons from the canonical transcript are indicated in grey in the plot below, 
            and labeled with their rank in that transcript.          
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
      initialPlotParametersSet: false
    }
  }, 
  methods: {
    isEven (number) {
      return number % 2 == 0
    },
    validateParameters () {
      this.errors = []

      if (!this.plotParameters.region) {
        this.errors.push('Region required.')
      }
      if (!this.plotParameters.windowSize) {
        this.errors.push('Window size required.')
      }
      if (!this.plotParameters.windowStride) {
        this.errors.push('Window stride required.')
      }   

      if (this.isEven(this.plotParameters.windowSize)) {
        this.errors.push('Window size must be odd.')
      }

      if (this.errors.length > 0) {
        return false
      } else { 
        return true
      }
    },
    getAPIData () {
      if ( this.validateParameters() ) {
        this.$store.dispatch('getMutationCounts', this.plotParameters)
        this.$store.dispatch('getCanonicalData', this.plotParameters.region)
      } else { 
        this.showSnackbar = true
      }
    },
  },
  computed: {
    ...mapState([
      'canonicalTranscript',
      'canonicalExons',
      'canonicalDataSet',
    ]),
    ...mapGetters([
      'fetchingAPIData'
    ]),
    canonicalTranscriptEnsembleUI () {
      return `https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=${this.canonicalTranscript.id}`
    },
    canonicalTranscriptExists () {
      return this.canonicalDataSet && Object.keys(this.canonicalTranscript).length > 0
    }
  },
  async created () { 
    this.plotParameters = await api.getInitialPlotParameters()
    this.initialPlotParametersSet = true
  }  
}
</script>

<style scoped> 
</style>

