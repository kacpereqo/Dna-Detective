<template>
  <div class="wrapper">
    <RnaInput @submit="submit" />
  </div>
</template>

<script>
import RnaInput from '@/components/RnaInput.vue';

export default {
  name: 'HomeView',
  components: {
    RnaInput
  }
}
</script>
