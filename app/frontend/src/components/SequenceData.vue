<template>
    <div class="sequence-data">
        <Visualization />
        <Propeties />
        <Charge />
        <!-- <Hydro /> -->
    </div>
</template>

<script>
import Visualization from '@/components/analize/Visualization.vue'
import Propeties from '@/components/analize/Propeties.vue'
import Hydro from '@/components/analize/Hydro.vue'
import Charge from '@/components/Charge.vue'

export default {
    name: 'SequenceData',
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
.sequence-data {
    padding: 1.5rem;
}

.property-wrapper {
    width: 100%;
    display: flex;
    flex-direction: column;
    text-align: center;
    justify-content: center;
    align-items: stretch;
}
</style>