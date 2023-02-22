import { createStore } from 'vuex'

export default createStore({
  state: {
    theme: localStorage.getItem('theme'),
    sequence: "",
  },
  getters: {
    theme: state => state.theme,
    sequence: state => state.sequence,
  },
  mutations: {
    toogleTheme(state) {
      if (state.theme === 'light-theme') {
        state.theme = 'dark-theme'
      }
      else {
        state.theme = 'light-theme'
      }
    },
    setSequence(state, sequence) {
      state.sequence = sequence;
    }
  },
  actions: {

  },
  modules: {
  }
})
