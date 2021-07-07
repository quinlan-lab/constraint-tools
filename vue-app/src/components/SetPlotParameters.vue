<!-- Based upon: https://vuematerial.io/components/form -->
<!-- sliders would be better for some inputs here, but the slider component is not yet implemented in this library -->
<!-- https://github.com/vuematerial/vue-material/blob/dev/ROADMAP.md -->
<!-- TODO: add form validation as per the website above -->
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

  </md-card>
</template> 

<script> 
import config from '@/assets/config'
import { mapState } from 'vuex'

export default {
  name: 'SetPlotParameters',
  data () {
    return {
      config: config,
    }
  }, 
  methods: {
    fetchData () {
      this.$store.dispatch('getMutationCounts', this.config)
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