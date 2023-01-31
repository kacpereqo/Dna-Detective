<template>
  <NavBar />
  <router-view />
</template>

<script>
import NavBar from '@/components/NavBar.vue'

export default {
  name: 'App',
  components: {
    NavBar
  }
}

</script>

<style>
@import url('https://fonts.googleapis.com/css2?family=Nunito:ital,wght@0,200;0,400;0,500;1,300&family=Roboto+Mono&display=swap');

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
}


body {
  min-height: 100vh;
  font-family: 'Nunito', sans-serif;
  box-sizing: border-box;
  padding: 0;
  margin: 0;
  overflow: auto;
}
</style>
