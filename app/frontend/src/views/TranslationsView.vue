
<template>
    <div class="wrapper">
        <!-- <Translations v-show="loaded" @loaded="loaded = true" /> -->
        <LoadingPrompt v-show="!loaded" :text="'Translacja'" :display_facts="true" />
    </div>
</template>

<script>
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
}
</script>

<style scoped>

</style>