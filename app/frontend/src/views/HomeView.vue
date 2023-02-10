<template>
  <div class="wrapper">
    <RnaInput @submit="submit" />
  </div>
</template>

<script>
import { useMeta } from 'vue-meta'
import RnaInput from '@/components/RnaInput.vue';

export default {
  name: 'HomeView',
  components: {
    RnaInput
  },
  setup() {
    useMeta({
      title: 'Strona główna',
    })
  },
}
</script>

<style scoped>
.wrapper {
  margin-top: -54px;
}
</style>
