
<template>
    <div class="wrapper">
        <Suspense>
            <Translations />
            <template #fallback>
                <LoadingPrompt :display_facts="true" />
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
        // clearIntervals() {
        // console.log('clear intervals');
        // }
    },
    setup() {
        useMeta({
            title: 'Translacja',
        })
    },
}
</script>

<style scoped></style>