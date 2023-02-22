import { createStore } from 'vuex'

export default createStore({
  state: {
    theme: localStorage.getItem('theme'),
    frame: localStorage.getItem('frame') || '',
  },
  getters: {
    theme: state => state.theme,
    frame: state => state.frame,
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
    setFrame(state, frame) {
      state.frame = frame;
      localStorage.setItem('frame', frame);
    }
  },
  actions: {

  },
  modules: {
  }
})
