<template>
    <div class="wrapper">
        <Visualization />
    </div>
</template>

<script>
import Visualization from '@/components/analize/Visualization.vue'

export default {
    name: 'AnalizeView',
    components: {
        Visualization
    },
    created() {
        this.id = this.$route.params.id;
    },
}
</script>

<style scoped>

</style>