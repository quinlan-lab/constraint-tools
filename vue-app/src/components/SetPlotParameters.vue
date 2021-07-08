<!-- Based upon: https://vuematerial.io/components/form -->
<!-- sliders would be better for some inputs here, but the slider component is not yet implemented in this library -->
<!-- https://github.com/vuematerial/vue-material/blob/dev/ROADMAP.md -->

<template> 
  <md-card 
    style="margin: 10px auto; " 
    class="md-layout-item md-size-75 md-xsmall-size-100"
  >
    <md-card-content>
      <div class="md-layout md-gutter">
        <div class="md-layout-item md-small-size-100">
          <md-field>
            <label for="region">Region</label>
            <md-input id="region" v-model="config.region" :disabled="fetchingData" />
          </md-field>
        </div>    

        <div class="md-layout-item md-small-size-100">
          <md-field>
            <label for="window-size">Window size</label>
            <md-input id="window-size" v-model="config.windowSize" :disabled="fetchingData" />
          </md-field>
        </div>    

        <div class="md-layout-item md-small-size-100">
          <md-field>
            <label for="window-stride">Window stride</label>
            <md-input id="window-stride" v-model="config.windowStride" :disabled="fetchingData" />
          </md-field>
        </div> 

      </div>
    </md-card-content>

    <md-card-actions>
      <md-button 
        v-on:click="fetchData" 
        class="md-icon-button md-primary"
        :disabled="fetchingData"
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
</template> 

<script> 
import config from '@/assets/config'
import { mapState } from 'vuex'

export default {
  name: 'SetPlotParameters',
  data () {
    return {
      errors: [],
      showSnackbar: false,
      config: config,
    }
  }, 
  methods:{
    isEven (number) {
      return number % 2 == 0
    },
    validateParameters () {
      this.errors = []

      if (!this.config.region) {
        this.errors.push('Region required.')
      }
      if (!this.config.windowSize) {
        this.errors.push('Window size required.')
      }
      if (!this.config.windowStride) {
        this.errors.push('Window stride required.')
      }
      
      if (this.isEven(this.config.windowSize)) {
        this.errors.push('Window size must be odd.')
      }

      if (this.errors.length > 0) {
        return false
      } else { 
        return true
      }
    },
    fetchData () {
      if ( this.validateParameters() ) {
        this.$store.dispatch('getMutationCounts', this.config)
      } else { 
        this.showSnackbar = true
      }
    },
  },
  computed: {
    ...mapState([
      'fetchingData',
    ])
  }
}
</script>

<style scoped> 
</style>