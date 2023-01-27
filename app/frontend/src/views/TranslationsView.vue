
<template>
    <div class="wrapper">
        <Translations />
    </div>
</template>

<script>
import Translations from '@/components/analize/Translation.vue';
import axios from 'axios';

export default {
    name: 'AnalizeView',
    components: {
        Translations
    },

}
</script>

<style scoped>

</style>