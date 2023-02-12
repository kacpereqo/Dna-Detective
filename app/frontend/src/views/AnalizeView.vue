<template>
    <div class="wrapper">
        <Sidebar @changeComponent="changeComponent" />
        <SequenceData :component="component" />
    </div>
</template>

<script>
import { useMeta } from 'vue-meta'

import Sidebar from '@/components/navigation/Sidebar.vue'
import SequenceData from '@/components/analize/SequenceData.vue'


export default {
    name: 'AnalizeView',
    components: {
        SequenceData,
        Sidebar,
    },
    data() {
        return {
            component: 'hydro',
        }
    },
    setup() {
        useMeta({
            title: 'Analiza',
        })
    },
    methods: {
        changeComponent(value) {
            this.component = value;
        },
    },
}
</script>

<style scoped>
.wrapper {
    display: flex;
    flex-direction: row;
    justify-content: flex-start;
    align-items: stretch;
}
</style>
