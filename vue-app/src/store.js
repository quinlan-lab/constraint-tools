import Vue from 'vue'
import Vuex from 'vuex'

// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Statements/import
import * as api from '@/api'

Vue.use(Vuex)

export const store = new Vuex.Store({
  state: {
    mutationCounts: null,
    fetchingMutationCounts: false,

    canonicalTranscript: null,
    canonicalExons: null,

    fetchingCanonicalData: false, 
    canonicalDataSet: false,
  },
  getters: {
    fetchingAPIData: state => {
      return state.fetchingMutationCounts || state.fetchingCanonicalData 
    }
  },
  mutations: {
    setMutationCounts (state, mutationCounts) {
      state.mutationCounts = mutationCounts
    },
    setFetchingMutationCounts (state, fetchingMutationCounts) {
      state.fetchingMutationCounts = fetchingMutationCounts
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
    async getMutationCounts ({ commit }, plotParameters) { 
      commit('setFetchingMutationCounts', true)
      commit('setMutationCounts', await api.getMutationCounts(plotParameters))
      commit('setFetchingMutationCounts', false)
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