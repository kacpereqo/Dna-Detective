<template>
    <div class="sequence-data">

        <Visualization v-show="component == 'visualization'" />
        <Propeties v-show="component == 'propeties'" />
        <Charge v-show="component == 'charge'" />
        <Hydro v-show="component == 'hydro'" />
    </div>
</template>

<script>
import SequenceName from '@/components/navigation/SequenceName.vue'
import Visualization from '@/components/analize/Visualization.vue'
import Propeties from '@/components/analize/Propeties.vue'
import Hydro from '@/components/analize/Hydro.vue'
import Charge from '@/components/analize/Charge.vue'

export default {
    name: 'SequenceData',
    components: {
        SequenceName,
        Visualization,
        Propeties,
        Hydro,
        Charge,
    },
    props: {
        component: {
            type: String,
            required: true
        }
    },
    data() {
        return {
            visibilityHistory: {
                visualization: false,
                propeties: false,
                charge: true,
                hydro: false,
            }
        }
    },
    watch: {
        component: function (value) {
            this.visibilityHistory[value] = true;
        }
    },
}

</script>

<style scoped>
.sequence-data {
    flex: 1;
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