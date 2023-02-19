<template>
    <select v-model="$i18n.locale" @change="saveLocale">
        <option v-for="(locale, i) in locales" :key="`locale-${i}`" :value="locale">
            {{ locale }}
        </option>
    </select>
</template>
 
<script>
export default {
    name: "LocaleSwitcher",
    data() {
        return { locales: ["pl", "en"] };
    },
    methods: {
        saveLocale() {
            localStorage.setItem('locale', this.$i18n.locale);
        }
    }
};
</script>