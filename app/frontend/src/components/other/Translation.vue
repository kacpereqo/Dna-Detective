<template>
    <div class="translations-wrapper">
        <h1>translacje</h1>
        <div class="translations">
            <ul>
                <li v-for="(translation, index) in translations">
                    <p>{{ translation.direction[0] }} {{ translation.frame }} {{ translation.direction[1] }}</p>
                    <p v-html="translation.openReadingFrames"></p>
                </li>
            </ul>
        </div>
    </div>
</template>

<script>
import axios from 'axios'
import { onMounted } from 'vue'
import { ref } from 'vue'
import { useRoute } from 'vue-router';

const getTranslation = async () => {
    const route = useRoute();
    const response = await axios.get(`http://127.0.0.1:8000/api/${route.params.id}/translate?is_reversed=true&is_forward=true`);
    return response.data;
}

export default {
    name: 'Translations',
    data() {
        return {
            id: '',
            translations: [],
        }
    },
    async setup() {
        const translations = await getTranslation();

        return {
            translations,
        }
    },
    methods: {
        getTranslation() {
            axios.get(`http://127.0.0.1:8000/api/7/translate?is_reversed=true&is_forward=true`)
                .then(res => {
                    this.translations = res.data;
                    this.$emit('loaded', true);
                })
        },

    }
}
</script>

<style>
a {
    color: var(--text-color);
}
</style>