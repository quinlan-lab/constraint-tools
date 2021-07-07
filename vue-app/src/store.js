import Vue from 'vue'
import Vuex from 'vuex'

import { api } from '@/api'

Vue.use(Vuex)

export const store = new Vuex.Store({
  state: {
    mutationCounts: null,
    fetchingData: false,
  },
  mutations: {
    setMutationCounts (state, mutationCounts) {
      state.mutationCounts = mutationCounts
    },
    setFetchingData (state, fetchingData) {
      state.fetchingData = fetchingData
    }
  },
  actions: {
    getMutationCounts ({ commit }, config) { 
      commit('setFetchingData', true)
      api(config).then(response => {
        commit('setMutationCounts', response.data)
        commit('setFetchingData', false)
      }).catch(error => {
        alert('getMutationCounts action: ' + error)
      })
    }
  }
})