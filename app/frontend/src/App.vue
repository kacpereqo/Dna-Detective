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
  display: flex;
  flex-direction: column;
  min-height: 100vh;
  flex: 1;
  justify-content: center;
  align-items: center;
}

#app {
  display: flex;
  flex-direction: column;
}

body {
  font-family: 'Nunito', sans-serif;
  box-sizing: border-box;
  padding: 0;
  margin: 0;
}
</style>
