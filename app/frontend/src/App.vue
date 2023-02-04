<template>
  <metainfo>
    <template v-slot:title="{ content }">{{ content? `${content} - DNA detective` : `SITE_NAME` }}</template>
  </metainfo>
  <NavBar />
  <router-view />
</template>

<script>
import { useMeta } from 'vue-meta'
import NavBar from '@/components/NavBar.vue'

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
  }

}

</script>

<style>
@import url('https://fonts.googleapis.com/css2?family=Nunito:ital,wght@0,200;0,400;0,500;1,300&family=Roboto+Mono&family=Roboto:wght@100;300;400;500;700&display=swap');

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
}


body {
  display: flex;
  flex-direction: column;
  min-height: 100vh;
  font-family: 'Nunito', sans-serif;
  box-sizing: border-box;
  padding: 0;
  margin: 0;

}
</style>
