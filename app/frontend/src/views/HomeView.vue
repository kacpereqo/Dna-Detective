<template>
  <div class="wrapper">
    <RnaInput @submit="submit" />
  </div>
</template>

<script>
import { useMeta } from 'vue-meta'
import RnaInput from '@/components/rnainput/RnaInput.vue';

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

</style>
