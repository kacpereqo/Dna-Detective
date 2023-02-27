import { createStore } from 'vuex'

export default createStore({
  state: {
    theme: localStorage.getItem('theme') || 'light-theme',
    frame: localStorage.getItem('frame') || '',
    jwt: localStorage.getItem('jwt') || '',
    isLogged: false,
    email: localStorage.getItem('email') || '',
    fontSize: localStorage.getItem('fontSize') || '16px',
  },
  getters: {
    theme: state => state.theme,
    frame: state => state.frame,
    jwt: state => state.jwt,
    isLogged: state => state.isLogged,
    email: state => state.email,
    fontSize: state => state.fontSize,
    sequence: state => state.sequence,
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
      localStorage.removeItem('email');
      state.isLogged = false;
    },
    setEmail(state, email) {
      state.email = email;
      localStorage.setItem('email', email);
    },
    setFontSize(state, fontSize) {
      document.documentElement.style.fontSize = fontSize + 'px';
      state.fontSize = fontSize;
      localStorage.setItem('fontSize', fontSize);
    },
  },
  actions: {

  },
  modules: {
  }
})
