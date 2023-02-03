
<template>
    <div class="wrapper">
        <Translations v-show="loaded" @loaded="loaded = true" />

        <span v-if="!loaded">
            <LoadingPrompt :text="'Translacja'" :display_facts="true" />
        </span>
    </div>
</template>

<script>
import { useMeta } from 'vue-meta'

import Translations from '@/components/analize/Translation.vue';
import LoadingPrompt from '@/components/analize/LoadingPrompt.vue';

export default {
    name: 'AnalizeView',
    components: {
        Translations,
        LoadingPrompt,
    },
    data() {
        return {
            loaded: false,
        }
    },
    methods: {
        destroy() {
            this.$destroy();
        }
    },
    setup() {
        useMeta({
            title: 'Translacja',
        })
    },
}
</script>

<style scoped>

</style>