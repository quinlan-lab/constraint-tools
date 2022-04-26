import Vue from 'vue'
import Vuex from 'vuex'

// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Statements/import
import * as api from '@/api'

Vue.use(Vuex)

export const store = new Vuex.Store({
  state: {
    expectedObservedCounts: null,
    fetchingExpectedObservedCounts: false,

    canonicalTranscript: null,
    canonicalExons: null,

    fetchingCanonicalData: false, 
    canonicalDataSet: false,
  },
  getters: {
    fetchingAPIData: state => {
      return state.fetchingExpectedObservedCounts || state.fetchingCanonicalData 
    }
  },
  mutations: {
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
    async getExpectedObservedCounts ({ commit }, plotParameters) { 
      commit('setFetchingExpectedObservedCounts', true)
      commit('setExpectedObservedCounts', await api.getExpectedObservedCounts(plotParameters))
      commit('setFetchingExpectedObservedCounts', false)
    },
    async getCanonicalData ({ commit }, region) { 
      // https://vuex.vuejs.org/guide/actions.html#composing-actions

      commit('setFetchingCanonicalData', true)

      const exons = await api.getExons(region)
      console.log('exons:')
      console.log(exons)      

      const transcriptIDs = [...new Set(exons.map(exon => exon.Parent))]
      console.log('transcriptIDs:')
      console.log(transcriptIDs)
      const canonicalTranscript = await api.getCanonicalTranscript(transcriptIDs)
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