<template>
    <select v-model="$i18n.locale" @change="saveLocale">
        <option v-for="(locale, i) in locales" :key="`locale-${i}`" :value="locale">
            {{ locale }}
        </option>
    </select>
</template>
 
<script>
import axios from 'axios';
export default {
    name: "LocaleSwitcher",
    data() {
        return { locales: ["pl", "en"] };
    },
    methods: {
        saveLocale() {
            localStorage.setItem('locale', this.$i18n.locale);
            this.updateUserPreference();
        },
        updateUserPreference() {
            if (this.$store.state.isLogged) {
                axios.post(`http://127.0.0.1:8000/api/preferences?key=lang&value=${this.$i18n.locale}`, {
                })
            }
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