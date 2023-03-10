<template>
  <metainfo>
    <template v-slot:title="{ content }">{{ content ? `${content} - DNA detective` : `SITE_NAME` }}</template>
  </metainfo>
  <NavBar />
  <keep-alive>
    <router-view />
  </keep-alive>
</template>

<script>
import { useMeta } from 'vue-meta'
import axios from 'axios'
import NavBar from '@/components/navigation/NavBar.vue'

export default {
  name: 'App',
  components: {
    NavBar
  },

  setup() {
    useMeta({
      title: 'NONE',
      htmlAttrs: { lang: 'pl-PL', amp: true }
    })
  },

  mounted() {
    const locale = localStorage.getItem('locale')
    if (locale) {
      this.$i18n.locale = locale
    }

    const theme = localStorage.getItem('theme')
    if (theme) {
      document.documentElement.className = theme
    }
    else {
      localStorage.setItem('theme', 'light-theme')
    }

    const fontSize = localStorage.getItem('font-size')
    if (fontSize) {
      document.documentElement.style.fontSize = fontSize
    }

    const jwt = localStorage.getItem('jwt')
    if (jwt) {
      this.$store.commit('setUser', jwt)
      axios.defaults.headers.common['Authorization'] = `Bearer ${jwt}`
    }

  }
}

</script>

<style>
@import url('https://fonts.googleapis.com/css2?family=Nunito:ital,wght@0,200;0,400;0,500;1,300&family=Roboto+Mono&family=Roboto:wght@100;300;400;500;700&display=swap');

:root {
  --visualization-filter: 0;
  --heading-color: #ffffff;
  --text-color: #000000;
  --accent-color-dark: #00000080;
  --accent-color: #00000040;
  --accent-color-light: #0000001a;
  --main-color: #5f5cff;
  --background-color: #ffffff;
  --chart-color: #0000ff4d;
  --main-color-opacity: #0400f5;
  --icon-filter: invert(0);
}

:root.dark-theme {
  --visualization-filter: invert(1) contrast(80%);
  --heading-color: #ffffff;
  --text-color: #ffffff;
  --accent-color-dark: #ffffff80;
  --accent-color: #ffffff40;
  --accent-color-light: #ffffff1a;
  --main-color: #7370f5b3;
  --main-color-opacity: #1511d4;
  --background-color: #191919;
  --icon-filter: invert(1);
  --chart-color: #6183ff33;
}

:root.high-contrast-theme {
  --visualization-filter: 0;
  --heading-color: #ffffff;
  --text-color: #000000;
  --accent-color-dark: #00000080;
  --accent-color: #00000040;
  --accent-color-light: #0000001a;
  --main-color: #5f5cff;
  --background-color: #ffffff;
  --chart-color: #0000ff4d;
  --main-color-opacity: #0400f5;
  --icon-filter: invert(0);
  --contrast: invert(1);
}

.wrapper {
  position: relative;
  display: flex;
  flex-direction: column;
  min-height: calc(100vh - 54px);
  height: fit-content;
  flex: 1;
  justify-content: center;
  align-items: center;
}

#app {
  position: relative;
  display: flex;
  flex-direction: column;
  min-height: 100vh;
  height: fit-content;
  overflow: hidden;
  background-color: var(--background-color);

}


body {
  color: var(--text-color);
  display: flex;
  overflow: auto;
  flex-direction: column;
  min-height: 100vh;
  font-family: 'Nunito', sans-serif;
  box-sizing: border-box;
  padding: 0;
  margin: 0;
}

a,
input,
textarea {
  color: var(--text-color);
}

* {
  transition: var(--transition);
}

html {
  filter: invert(1);
  filter: var(--contrast);
  font-size: 16px;
}
</style>
