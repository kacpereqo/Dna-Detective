
<template>
    <div class="wrapper">
        <Suspense :display_facts="true">
            <Translations />
            <template #fallback>
                <LoadingPrompt />
            </template>
        </Suspense>
    </div>
</template>


<script>
import { useMeta } from 'vue-meta'

import Translations from '@/components/other/Translation.vue';
import LoadingPrompt from '@/components/other/LoadingPrompt.vue';

export default {
    name: 'AnalizeView',
    components: {
        LoadingPrompt,
        Translations,
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