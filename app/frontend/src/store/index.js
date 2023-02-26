import { createStore } from 'vuex'

export default createStore({
  state: {
    theme: localStorage.getItem('theme'),
    frame: localStorage.getItem('frame') || '',
    jwt: localStorage.getItem('jwt') || '',
  },
  getters: {
    theme: state => state.theme,
    frame: state => state.frame,
    jwt: state => state.jwt,
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
    },
    setUser(state, jwt) {
      state.jwt = jwt;
    }
  },
  actions: {

  },
  modules: {
  }
})
