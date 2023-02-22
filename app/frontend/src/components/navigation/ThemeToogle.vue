<template>
    <div class="theme-switch">
        <img src="@/assets/contrast.svg" @click="toogleTheme">
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

            document.documentElement.style.setProperty('--transition', '0.5s');
            setTimeout(() => {
                document.documentElement.style.setProperty('--transition', '0s');
            }, 400);

            if (this.theme === 'light-theme') {
                this.theme = 'dark-theme';
            } else if (this.theme === 'dark-theme') {
                this.theme = 'high-contrast-theme';
            } else if (this.theme === 'high-contrast-theme') {
                this.theme = 'light-theme';
            }


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

<style scoped lang="scss">
.theme-switch {
    display: flex;
    align-items: center;
    justify-content: center;
    margin: 0 1rem;
    cursor: pointer;
    filter: var(--icon-filter);
}
</style>