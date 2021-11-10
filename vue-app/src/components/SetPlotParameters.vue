<!-- Based upon: https://vuematerial.io/components/form -->
<!-- sliders would be better for some inputs here, but the slider component is not yet implemented in this library -->
<!-- https://github.com/vuematerial/vue-material/blob/dev/ROADMAP.md -->

<template> 
  <div class="md-layout md-gutter">
    <div class="md-layout-item md-small-size-100">
      <md-progress-bar v-if="!initialPlotParametersSet" md-mode="indeterminate" />
      <md-card 
        v-else
        style="margin: 10px auto;" 
        class="md-layout-item md-size-75 md-xsmall-size-100"
      >

        <md-card-content>
          <div class="md-layout md-gutter">
            <div class="md-layout-item md-small-size-100">
              <md-field>
                <label for="region">Region</label>
                <md-input id="region" v-model="plotParameters.region" :disabled="fetchingAPIData" />
              </md-field>
            </div>    

            <div class="md-layout-item md-small-size-100">
              <md-field>
                <label for="window-size">Window size</label>
                <md-input id="window-size" v-model="plotParameters.windowSize" :disabled="fetchingAPIData" />
              </md-field>
            </div>    

            <div class="md-layout-item md-small-size-100">
              <md-field>
                <label for="window-stride">Window stride</label>
                <md-input id="window-stride" v-model="plotParameters.windowStride" :disabled="fetchingAPIData" />
              </md-field>
            </div> 
          </div>
        </md-card-content>

        <md-card-actions>
          <md-button 
            v-on:click="getAPIData" 
            class="md-icon-button md-primary"
            :disabled="fetchingAPIData"
          >
            <md-icon>refresh</md-icon>
          </md-button>      
        </md-card-actions>

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

    <div 
      v-if="canonicalDataSet && Object.keys(canonicalTranscript).length !== 0"
      class="md-layout-item md-small-size-100"
    >
      <md-card 
        style="margin: 10px auto;" 
        class="md-layout-item md-size-75 md-xsmall-size-100"
      >
        <md-card-content>
          <!-- https://www.creative-tim.com/vuematerial/components/table -->
          <md-table v-if="canonicalDataSet">
            <md-table-row>
              <md-table-head>Canonical Transcript Property</md-table-head>
              <md-table-head>Value</md-table-head>
            </md-table-row>

            <md-table-row href="https://www.google.com" target="_blank">
              <md-table-cell>Strand</md-table-cell>
              <md-table-cell md-numeric>{{ canonicalTranscript.strand }}</md-table-cell>
            </md-table-row>

            <md-table-row>
              <md-table-cell>Name</md-table-cell>
              <md-table-cell>{{ canonicalTranscript.display_name }}</md-table-cell>
            </md-table-row>

            <md-table-row>
              <md-table-cell>Ensemble ID</md-table-cell>
              <md-table-cell>
                <a :href="canonicalTranscriptEnsembleUI" target="_blank">
                  {{ canonicalTranscript.id }}
                </a>                        
              </md-table-cell>
            </md-table-row>
          </md-table>  
          Exons from the canonical transcript are indicated in grey in the plot below, 
          and labeled with their rank in that transcript.
        </md-card-content> 
      </md-card>
    </div> 
    
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

