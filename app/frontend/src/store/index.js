import { createStore } from 'vuex'

export default createStore({
  state: {
    theme: localStorage.getItem('theme'),
    frame: localStorage.getItem('frame') || '',
    jwt: localStorage.getItem('jwt') || '',
    isLogged: false,
  },
  getters: {
    theme: state => state.theme,
    frame: state => state.frame,
    jwt: state => state.jwt,
    isLogged: state => state.isLogged,
  },
  mutations: {
    toogleTheme(state, theme) {
      state.theme = theme;
      document.documentElement.className = theme;
      localStorage.setItem('theme', theme);
    },
    setFrame(state, frame) {
      state.frame = frame;
      localStorage.setItem('frame', frame);
    },
    setUser(state, jwt) {
      state.jwt = jwt;
      localStorage.setItem('jwt', jwt);
      state.isLogged = true;
    },
    logout(state) {
      state.jwt = '';
      localStorage.removeItem('jwt');
      state.isLogged = false;
    }
  },
  actions: {

  },
  modules: {
  }
})
