import Vue from 'vue'
import Vuex from 'vuex'

// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Statements/import
import * as api from '@/api'

Vue.use(Vuex)

export const store = new Vuex.Store({
  state: {
    selectedGenomicPosition: null, 
    distributionN: null,
    fetchingDistributionN: false,

    modelParameters: null, 
    fetchingModelParameters: false, 
    modelParametersSet: false,
    
    neutralRegions: null, 
    fetchingNeutralRegions: false, 
    neutralRegionsSet: false,

    expectedObservedCounts: null,
    fetchingExpectedObservedCounts: false,

    canonicalTranscript: null,
    canonicalExons: null,

    fetchingCanonicalData: false, 
    canonicalDataSet: false,
  },
  getters: {
    fetchingTimeSeriesData: state => {
      return (
        state.fetchingExpectedObservedCounts || 
        state.fetchingCanonicalData || 
        state.fetchingNeutralRegions 
      )
    },
    fetchingDistributionData: state => {
      return (
        state.fetchingDistributionN // || state.fetchingDistributionK
      )
    }
  },
  mutations: {
    setSelectedGenomicPosition (state, selectedGenomicPosition) {
      state.selectedGenomicPosition = selectedGenomicPosition
    },
    setDistributionN (state, distributionN) {
      state.distributionN = distributionN
    },
    setFetchingDistributionN (state, fetchingDistributionN) {
      state.fetchingDistributionN = fetchingDistributionN
    },

    setModelParameters (state, modelParameters) {
      state.modelParameters = modelParameters
      state.modelParametersSet = true
    },
    setFetchingModelParameters (state, fetchingModelParameters) {
      state.fetchingModelParameters = fetchingModelParameters
    },

    setNeutralRegions (state, neutralRegions) {
      state.neutralRegions = neutralRegions
      state.neutralRegionsSet = true
    },
    setFetchingNeutralRegions (state, fetchingNeutralRegions) {
      state.fetchingNeutralRegions = fetchingNeutralRegions
    },

    setExpectedObservedCounts (state, expectedObservedCounts) {
      state.expectedObservedCounts = expectedObservedCounts
    },
    setFetchingExpectedObservedCounts (state, fetchingExpectedObservedCounts) {
      state.fetchingExpectedObservedCounts = fetchingExpectedObservedCounts
    },

    setCanonicalTranscript (state, canonicalTranscript) {
      state.canonicalTranscript = canonicalTranscript
    },
    setCanonicalExons (state, canonicalExons) {
      state.canonicalExons = canonicalExons
    },

    setFetchingCanonicalData (state, fetchingCanonicalData) {
      state.fetchingCanonicalData = fetchingCanonicalData
    },
    setCanonicalDataSet (state, canonicalDataSet) {
      state.canonicalDataSet = canonicalDataSet
    },
  },
  actions: {
    async getDistributionN ({ commit }, window) { 
      commit('setFetchingDistributionN', true)
      commit('setDistributionN', await api.getDistributionN(window))
      commit('setFetchingDistributionN', false)
    },
    async getModelParameters ({ commit }) { 
      commit('setFetchingModelParameters', true)
      commit('setModelParameters', await api.getModelParameters())
      commit('setFetchingModelParameters', false)
    },
    async getNeutralRegions ({ commit }, plotParameters) { 
      commit('setFetchingNeutralRegions', true)
      commit('setNeutralRegions', await api.getNeutralRegions(plotParameters))
      commit('setFetchingNeutralRegions', false)
    },
    async getExpectedObservedCounts ({ commit }, plotParameters) { 
      commit('setFetchingExpectedObservedCounts', true)
      commit('setExpectedObservedCounts', await api.getExpectedObservedCounts(plotParameters))
      commit('setFetchingExpectedObservedCounts', false)
    },
    async getCanonicalData ({ commit }, payload) { 
      // https://vuex.vuejs.org/guide/actions.html#composing-actions

      commit('setFetchingCanonicalData', true)

      const exons = await api.getExons(payload.region, payload.genomeBuild)
      console.log('exons:')
      console.log(exons)      

      const transcriptIDs = [...new Set(exons.map(exon => exon.Parent))]
      console.log('transcriptIDs:')
      console.log(transcriptIDs)
      const canonicalTranscript = await api.getCanonicalTranscript(transcriptIDs, payload.genomeBuild)
      console.log('canonicalTranscript:')
      console.log(canonicalTranscript)
      commit('setCanonicalTranscript', canonicalTranscript)

      const canonicalExons = exons.filter(exon => exon.Parent == canonicalTranscript.id)
      console.log('canonicalExons')
      console.log(canonicalExons)
      commit('setCanonicalExons', canonicalExons)

      commit('setFetchingCanonicalData', false)
      commit('setCanonicalDataSet', true)
    }
  }
})