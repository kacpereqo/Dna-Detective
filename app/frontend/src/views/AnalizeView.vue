<template>
    <div class="wrapper">
        <Sidebar />
        <div class="sequence-data">
            <Visualization />
            <Propeties />
            <Charge />
            <Hydro />
            <Hydro />
            <Hydro />
            <Hydro />
            <Hydro />
            <Hydro />
            <Hydro />
            <Hydro />
            <Hydro />
        </div>
    </div>
</template>

<script>
import Visualization from '@/components/analize/Visualization.vue'
import Propeties from '@/components/analize/Propeties.vue'
import Hydro from '@/components/analize/Hydro.vue'
import Charge from '@/components/analize/Charge.vue'

import Sidebar from '@/components/Sidebar.vue'


export default {
    name: 'AnalizeView',
    components: {
        Visualization,
        Propeties,
        Hydro,
        Charge,
        Sidebar,
    },
    created() {
        this.id = this.$route.params.id;
    },
}
</script>

<style scoped>
.wrapper {
    display: flex;
    flex-direction: row;
    align-items: stretch;
}

.sequence-data {
    display: flex;
    align-items: center;
    justify-content: center;
    flex-direction: column;
    flex: 1;
}
</style>
