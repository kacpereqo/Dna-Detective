<template>
    <div class="nav-wrapper">
        <nav>
            <div class="logo">
                <router-link to="/">
                </router-link>
            </div>
            <div class="nav-links">
                <div class='desktop'>
                    <ul>
                        <li style="flex:1;">
                            <div class="main-links">
                                <ul>
                                    <li><router-link to="/#" class="nav-link">Strona główna</router-link></li>
                                    <li><router-link to="/#" class="nav-link">Analiza</router-link></li>
                                    <li>
                                        <DropDownMenu :hrefs="[{ name: 'Tabela 1', path: '/' },
                                        { name: 'Tablea 2', path: '/' },]" :name="'Tabele'" />
                                    </li>
                                    <li class="drop-down">
                                        <DropDownMenu :hrefs="[{ name: 'Narzedzie 1', path: '/' },
                                        { name: 'Narzedzie 2', path: '/' },]" :name="'Narzędzia'" />
                                    </li>
                                </ul>
                            </div>
                        </li>
                        <li>
                            <div class="left">
                                <ul>
                                    <li class="toogle-switch">
                                        <ThemeToogle />
                                    </li>
                                    <li>
                                        <DropDownMenu :hrefs="[{ name: 'Zaloguj', path: '/' },
                                        { name: 'Zarejestruj', path: '/' }, { name: 'Preferencje' }]"
                                            :name="'Konto'" />
                                    </li>
                                </ul>
                            </div>
                        </li>
                    </ul>
                </div>
                <div class="mobile">
                    <div class="icon" @click="showMenu = !showMenu"
                        :style="{ '--width-top': showMenu ? '1rem' : '2rem', '--width-botom': showMenu ? '2rem' : '1rem' }">
                        <span>
                        </span>
                    </div>
                    <div class="menu" v-if="showMenu">
                        <ul>
                            <li>
                                <router-link to="/#"> Strona główna</router-link>
                            </li>

                        </ul>
                    </div>
                    <img src="@/assets/user.svg" class="mobile-user">
                </div>
            </div>
        </nav>
    </div>
</template>

<script>
import DropDownMenu from './DropDownMenu.vue';
import ThemeToogle from './ThemeToogle.vue';

export default {
    name: 'NavBar',
    data() {
        return {
            showMenu: false
        }
    },
    components: {
        DropDownMenu,
        ThemeToogle
    },

};
</script>
<style scoped>
nav {
    z-index: 2;
    position: relative;
    display: flex;
    flex-direction: row;
    line-height: 54px;
    box-shadow: 0 1px 1px var(--accent-color), 0 3px 2px var(--accent-color-light);
    height: 54px;
}

p {
    margin: 0;
}

a {
    color: var(--text-color);
    text-decoration: none;
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


/* mobile */

@media only screen and (max-width: 960px) {
    .desktop {
        display: none;
    }

    .mobile {
        display: flex !important;
    }
}

/* desktop  */

@media only screen and (max-width: 1264px) {

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

a:hover {
    background-color: var(--accent-color);
}

.nav-links ul {
    display: flex;
    justify-content: center;
    align-items: center;
}

.main-links {
    flex: 1;
}

.nav-links {
    flex: 1;
}

.left {
    margin-left: auto;
}


.mobile {
    height: 100%;
    z-index: 5;
    justify-content: flex-end;
    align-items: center;
    display: none;
}

.icon {
    position: relative;
    width: 2rem;
    height: 1rem;
    cursor: pointer;
}


.icon span {
    position: absolute;
    top: 50%;
    transform: translateY(-50%);
    right: 0;
    width: 1.5rem;
    border-top: 2px solid var(--text-color);
}

.icon::after {
    content: '';
    width: 2rem;
    position: absolute;
    bottom: 0;
    right: 0;
    border-top: 2px solid var(--text-color);
}

.icon::before {
    content: '';
    width: 1rem;
    position: absolute;
    top: 0;
    right: 0;
    border-top: 2px solid var(--text-color);
}


.icon {
    position: relative;
    width: 2rem;
    height: 1rem;
    cursor: pointer;
}


.icon span {
    position: absolute;
    top: 50%;
    transform: translateY(-50%);
    right: 0;
    width: 1.5rem;
    border-top: 2px solid var(--text-color);
}

.icon::after {
    content: '';
    width: var(--width-botom);
    position: absolute;
    bottom: 0;
    right: 0;
    border-top: 2px solid var(--text-color);
    transition: 0.2s;
}

.icon::before {
    content: '';
    width: var(--width-top);
    position: absolute;
    top: 0;
    right: 0;
    border-top: 2px solid var(--text-color);
    transition: 0.2s;
}

.menu {
    background-color: var(--background-color);
    position: absolute;
    width: 100vw;
    left: 0;
    height: fit-content;
    top: 100%;
}


.menu ul {
    width: 100%;
    height: 100%;
    margin: 0;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
}


.menu ul li {
    width: 100%;
    text-align: center;
    border: 1px solid var(--accent-color);
}

.mobile-user {
    width: 1.5rem;
    height: 1.5rem;
    margin: 0 0.5rem;
    padding: 0.5rem;
    cursor: pointer;
    filter: var(--icon-filter);
}
</style>