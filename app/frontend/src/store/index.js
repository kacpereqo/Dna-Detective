import { createStore } from 'vuex'

export default createStore({
  state: {
    theme: localStorage.getItem('theme'),
  },
  getters: {
    theme: state => state.theme,
  },
  mutations: {
    toogleTheme(state) {
      if (state.theme === 'light-theme') {
        state.theme = 'dark-theme'
      }
      else {
        state.theme = 'light-theme'
      }
    }
  },
  actions: {

  },
  modules: {
  }
})
