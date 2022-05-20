import Vue from 'vue'
import Vuex from 'vuex'

// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Statements/import
import * as api from '@/api'

Vue.use(Vuex)

export const store = new Vuex.Store({
  state: {
    observedValues: null,  
    selectedGenomicPosition: null, 
    distributions: null,
    fetchingDistributions: false,

    modelParameters: null, 
    fetchingModelParameters: false, 
    modelParametersSet: false,
    
    neutralRegions: null, 
    fetchingNeutralRegions: false, 
    neutralRegionsSet: false,

    sequenceData: null, 
    fetchingSequenceData: false, 
    sequenceDataSet: false,

    expectedObservedCounts: null,
    fetchingExpectedObservedCounts: false,

    canonicalTranscripts: null,
    canonicalExons: null,

    fetchingCanonicalData: false, 
    canonicalDataSet: false,

    exonColor: 'rgba(255, 0, 0, 0.3)',
    neutralRegionColor: 'rgba(0, 0, 0, 0.1)', // '#d3d3d3', 
  },
  getters: {
    fetchingTimeSeriesData: state => {
      return (
        state.fetchingExpectedObservedCounts || 
        state.fetchingCanonicalData || 
        state.fetchingNeutralRegions ||
        state.fetchingSequenceData 
      )
    },
  },
  mutations: {
    setObservedValues (state, observedValues) {
      state.observedValues = observedValues
    },
    setSelectedGenomicPosition (state, selectedGenomicPosition) {
      state.selectedGenomicPosition = selectedGenomicPosition
    },
    setDistributions (state, distributions) {
      state.distributions = distributions
    },
    setFetchingDistributions (state, fetchingDistributions) {
      state.fetchingDistributions = fetchingDistributions
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

    setSequenceData (state, sequenceData) {
      state.sequenceData = sequenceData
      state.sequenceDataSet = true
    },
    setFetchingSequenceData (state, fetchingSequenceData) {
      state.fetchingSequenceData = fetchingSequenceData
    },

    setExpectedObservedCounts (state, expectedObservedCounts) {
      state.expectedObservedCounts = expectedObservedCounts
    },
    setFetchingExpectedObservedCounts (state, fetchingExpectedObservedCounts) {
      state.fetchingExpectedObservedCounts = fetchingExpectedObservedCounts
    },

    setCanonicalTranscripts (state, canonicalTranscripts) {
      state.canonicalTranscripts = canonicalTranscripts
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
    async getDistributions ({ commit }, payload) { 
      commit('setFetchingDistributions', true)
      commit('setDistributions', await api.getDistributions(payload))
      commit('setFetchingDistributions', false)
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
    async getSequenceData ({ commit }, plotParameters) { 
      commit('setFetchingSequenceData', true)
      commit('setSequenceData', await api.getSequenceData(plotParameters))
      commit('setFetchingSequenceData', false)
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
      const canonicalTranscripts = await api.getCanonicalTranscripts(transcriptIDs, payload.genomeBuild)
      console.log('canonicalTranscripts:')
      console.log(canonicalTranscripts)
      commit('setCanonicalTranscripts', canonicalTranscripts)

      const canonicalTranscriptIDs = new Set(canonicalTranscripts.map(canonicalTranscript => canonicalTranscript.id))
      const canonicalExons = exons.filter(exon => canonicalTranscriptIDs.has(exon.Parent))
      console.log('canonicalExons')
      console.log(canonicalExons)
      commit('setCanonicalExons', canonicalExons)

      commit('setFetchingCanonicalData', false)
      commit('setCanonicalDataSet', true)
    }
  }
})