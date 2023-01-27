<template>
    <div class="wrapper">
        <Charge />
        <Hydro />
        <Propeties />
        <Visualization />
    </div>
</template>

<script>
import Visualization from '@/components/analize/Visualization.vue'
import Propeties from '@/components/analize/Propeties.vue'
import Hydro from '@/components/analize/Hydro.vue'
import Charge from '@/components/analize/Charge.vue'


export default {
    name: 'AnalizeView',
    components: {
        Visualization,
        Propeties,
        Hydro,
        Charge,
    },
    created() {
        this.id = this.$route.params.id;
    },
}
</script>

<style scoped>

</style>