<template>
    <div class="nav-wrapper">
        <nav>
            <div class="logo">
                <router-link to="/">
                </router-link>
            </div>
            <div class="nav-links">
                <ul>
                    <li><router-link to="/#" class="nav-link">Strona główna</router-link></li>
                    <li><router-link to="/#" class="nav-link">Analiza</router-link></li>
                    <li>
                        <DropDownMenu :hrefs="[{ name: 'Tabela 1', href: '/#' },
                        { name: 'Tablea 2', href: '/#' },]" :name="'Tabele'" />
                    </li>
                    <li class="drop-down">
                        <DropDownMenu :hrefs="[{ name: 'Narzedzie 1', href: '/#' },
                        { name: 'Narzedzie 2', href: '/#' },]" :name="'Narzędzia'" />
                    </li>
                </ul>
            </div>
            <div class="login">
                <ul>
                    <li class="toogle-switch" @click="toggleDarkMode">Darkmode</li>
                    <li>
                        <DropDownMenu :hrefs="[{ name: 'Zaloguj', href: '/#' },
                        { name: 'Zarejestruj', href: '/#' },]" :name="'Konto'" />
                    </li>
                </ul>
            </div>
        </nav>
    </div>
</template>

<script>
import DropDownMenu from './DropDownMenu.vue';

export default {
    name: 'NavBar',
    components: {
        DropDownMenu
    },
    methods: {
        toggleDarkMode() {
            this.$store.commit('toggleDarkMode');
        }
    }
};
</script>
<style scoped>
.toogle-switch {
    cursor: pointer;
}

nav {
    z-index: 4;
    position: relative;
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    line-height: 54px;
    padding: 0 1rem;
    box-shadow: 0 2px 4px var(--accent-color);
    height: 54px;
}

p {
    margin: 0;
}

.logo a {
    display: block;
    height: 48px;
    background-image: url('../assets/logo.png');
    background-size: contain;
    background-repeat: no-repeat;
    background-position: center;
    aspect-ratio: 4/1;
}

@media screen and (max-width: 1100px) {
    .logo a {
        background-image: url('@/assets/small-logo.svg');
        background-size: contain;
        background-repeat: no-repeat;
        background-position: center;
        aspect-ratio: 1/1;
    }
}

.logo {
    position: relative;
}

ul {
    display: flex;
    justify-content: center;
    align-items: center;
    flex-direction: row;
    flex: 1;
    margin: 0;
    padding: 0;
}

ul li {
    list-style: none;
}

ul li {
    padding: 0 0.75rem;
}


a {
    padding: 0 0.4rem;
    text-decoration: none;
    display: block;
}

.login {
    display: flex;
    justify-content: center;
    align-items: center;
    position: relative;
}

a:hover {
    background-color: var(--accent-color);
}
</style>