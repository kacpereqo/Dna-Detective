import { createStore } from 'vuex'

export default createStore({
  state: {
    isDarkMode: false,
  },
  getters: {
    isDarkMode: state => state.isDarkMode,
  },
  mutations: {
    toggleDarkMode(state) {
      state.isDarkMode = !state.isDarkMode;
    }
  },
  actions: {

  },
  modules: {
  }
})
