<template>
    <div class="translations-wrapper">
        <h1>translacje</h1>
        <div class="translations">
            <ul>
                <li v-for="(translation, index) in translations">
                    <p>{{ translation.direction[0] }} {{ translation.frame }} {{ translation.direction[1] }}</p>
                    <!-- <div>{{ translation.translatedFrame }}</div> -->
                    <p v-html="translation.openReadingFrames"></p>
                </li>
            </ul>
        </div>
    </div>
</template>

<script>
import axios from 'axios'

export default {
    name: 'Translations',
    data() {
        return {
            hrefs: [],
            id: '',
            translations: [],
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getTranslation();
    },
    methods: {
        getTranslation() {
            axios.get(`http://127.0.0.1:8000/api/${this.id}/translate?is_reversed=true&is_forward=true`)
                .then(res => {
                    this.translations = res.data;
                    this.$emit('loaded', true);
                })
        },

    }
}
</script>