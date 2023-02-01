<template>
    <div class="wrapper">
        <Sidebar />
        <SequenceData />
    </div>
</template>

<script>

import Sidebar from '@/components/Sidebar.vue'
import SequenceData from '@/components/SequenceData.vue'


export default {
    name: 'AnalizeView',
    components: {
        SequenceData,
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
    justify-content: flex-start;
    align-items: stretch;
}
</style>
