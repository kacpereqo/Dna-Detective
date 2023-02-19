<template>
    <div class="toogle-wrapper">
        <input id="chck" type="checkbox" ref="chck" @click="toogleTheme">
        <label for="chck" class="check-trail">
            <span class="check-handler"></span>
        </label>
    </div>
</template>

<script>

export default {
    name: 'ThemeToogle',
    methods: {
        toogleTheme() {
            const theme = localStorage.getItem('theme');

            document.documentElement.style.setProperty('--transition', '0.4s');
            setTimeout(() => {
                document.documentElement.style.setProperty('--transition', '0s');
            }, 400);

            if (theme === 'dark-theme') {
                localStorage.setItem('theme', 'light-theme');
                document.documentElement.className = 'light-theme';
            } else {
                localStorage.setItem('theme', 'dark-theme');
                document.documentElement.className = 'dark-theme';
            }
            this.$store.commit('toogleTheme', theme);
        }
    },
    mounted() {
        const theme = localStorage.getItem('theme');
        if (theme === 'dark-theme') {
            this.$refs.chck.checked = true;
        } else {
            document.documentElement.className = 'light-theme';
        }
    }
}
</script>

<style scoped lang="scss">
.toogle-wrapper {
    cursor: pointer;
}

$transition: all .5s ease;

input[type="checkbox"] {
    position: absolute;
    opacity: 0;
    z-index: -1;
}

.check-trail {
    display: flex;
    align-items: center;
    width: 4rem;
    height: 2rem;
    background: transparent;
    border: 1px solid var(--accent-color-dark);
    border-radius: 2.5em;
    transition: $transition;
    cursor: pointer;
}

.check-handler {
    display: flex;
    margin-left: 0.25rem;
    justify-content: center;
    align-items: center;
    width: 1.5rem;
    height: 1.5rem;
    background: white;
    border-radius: 50%;
    transition: $transition;
    border: 1px solid var(--accent-color-dark);

    &:before {
        content: "";
        background-image: url('../../assets/sun.svg');
        background-size: contain;
        background-repeat: no-repeat;
        background-position: center;
        width: 1.5rem;
        height: 1.5rem;

    }
}

input[type="checkbox"]:checked+.check-trail {

    .check-handler {
        margin-left: 50%;
        background: transparent;

        &::before {
            content: "";
            background-image: url('../../assets/moon.svg');
            background-size: contain;
            background-repeat: no-repeat;
            background-position: center;
            width: 1.5rem;
            height: 1.5rem;
            filter: invert(1);
        }
    }
}
</style>