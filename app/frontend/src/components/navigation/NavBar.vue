<template>
    <div class="nav-wrapper">
        <nav>
            <div class="logo">
                <router-link to="/">
                    <div class="logo-img"></div>
                </router-link>
            </div>
            <div class="nav-links">
                <div class='desktop'>
                    <ul>
                        <li style="flex:1;">
                            <div class="main-links">
                                <ul>
                                    <li><router-link to="/#" class="nav-link">{{ $t('navbar.home') }}</router-link></li>
                                    <li><router-link to="/#" class="nav-link">{{ $t('navbar.analise') }}</router-link>
                                    </li>
                                    <li><router-link v-if="isLogged" to="/profile" class="nav-link">Historia</router-link>
                                    </li>

                                </ul>
                            </div>
                        </li>
                        <li>
                            <div class="left">
                                <ul>
                                    <li>
                                        <FontSwitch />
                                    </li>
                                    <li class="toogle-switch">
                                        <ThemeToogle />
                                    </li>
                                    <li>
                                        <LocaleSwitcher />
                                    </li>
                                    <li>
                                        <DropDownMenu :hrefs="isLogged ? [{ name: $t('navbar.account.logout'), path: '/logout' }] : [{ name: $t('navbar.account.login'), path: '/login' },
                                        { name: $t('navbar.account.register'), path: '/register' }]"
                                            :name="$t('navbar.account.account')" />
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
                                <router-link to="/#"> {{ $t('navbar.home') }}</router-link>
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
import LocaleSwitcher from './LocaleSwitcher.vue';
import FontSwitch from './FontSwitch.vue';
import { mapState } from 'vuex'

export default {
    name: 'NavBar',
    data() {
        return {
            showMenu: false,
        }
    },
    components: {
        LocaleSwitcher,
        DropDownMenu,
        ThemeToogle,
        FontSwitch
    },
    computed: {
        ...mapState({
            isLogged: 'isLogged'
        })

    },


};
</script>
<style scoped>
nav {
    z-index: 2;
    position: relative;
    display: flex;
    flex-direction: row;
    line-height: 3.375rem;
    box-shadow: 0 1px 1px var(--accent-color), 0 3px 2px var(--accent-color-light);
    height: 3.375rem;
}

p {
    margin: 0;
}

a {
    color: var(--text-color);
    text-decoration: none;
}

.logo a {
    padding: 0 1rem;
    height: 3.375rem;
    aspect-ratio: 4/1;
}

.logo-img {
    filter: var(--icon-filter);
    width: 100%;
    height: 100%;
    /* position: absolute; */
    background-image: url('../../assets/logo.png');
    background-size: contain;
    background-repeat: no-repeat;
    background-position: center;
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

    .logo-img {
        background-image: url('@/assets/small-logo.svg');
        width: 3.375rem;
        height: 3.375rem;
    }

    .logo a {
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
    padding: 0 0.15rem;
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
    margin-right: 1.75rem;
}

.left li {
    margin-right: 0.5rem;
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