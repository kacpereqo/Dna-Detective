<template>
    <div class="theme-switch">
        <select v-model="theme" @change="toogleTheme">
            <option value="light-theme">Light</option>
            <option value="dark-theme">Dark</option>
            <option value="high-contrast-theme">High Contrast</option>
        </select>
    </div>
</template>

<script>

export default {
    name: 'ThemeToogle',
    data() {
        return {
            theme: '',
        }
    },
    methods: {
        toogleTheme() {

            document.documentElement.style.setProperty('--transition', '0.4s');
            setTimeout(() => {
                document.documentElement.style.setProperty('--transition', '0s');
            }, 400);

            document.documentElement.className = this.theme;
            localStorage.setItem('theme', this.theme);
            this.$store.commit('toogleTheme', this.theme);
        }
    },
    mounted() {
        this.theme = localStorage.getItem('theme');

    }
}
</script>

<style scoped lang="scss"></style>