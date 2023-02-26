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

<style scoped>
select {
    background-color: var(--background-color);
    color: var(--text-color);
    border: 1px solid var(--accent-color-dark);
    border-radius: 5px;
    padding: 0.5rem;
    font-size: 0.9rem;
    cursor: pointer;
    outline: none;
}

select:hover {
    filter: brightness(0.85);
}
</style>