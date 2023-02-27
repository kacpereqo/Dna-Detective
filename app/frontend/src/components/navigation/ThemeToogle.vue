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
            let theme = localStorage.getItem('theme');
            document.documentElement.style.setProperty('--transition', '0.5s');
            setTimeout(() => {
                document.documentElement.style.setProperty('--transition', '0s');
            }, 400);

            if (theme === 'light-theme') {
                theme = 'dark-theme';
            } else if (theme === 'dark-theme') {
                theme = 'high-contrast-theme';
            } else if (theme === 'high-contrast-theme') {
                theme = 'light-theme';
            }

            this.$store.commit('toogleTheme', theme);
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
    cursor: pointer;
    filter: var(--icon-filter);
}
</style>