<template>
    <div class="nav-wrapper">
        <nav>
            <div class="logo">
                <router-link to="/">
                </router-link>
            </div>
            <div class="nav-links">
                <ul>
                    <li><router-link to="/" class="nav-link">Strona główna</router-link></li>
                    <li><router-link to="/#" class="nav-link">Analiza</router-link></li>
                    <li class="drop-down">
                        <p>Tabele</p>
                        <div class="drop-down-content">
                            <router-link to="/#">Tabela 1</router-link>
                            <router-link to="/#">Tabela 2</router-link>
                            <router-link to="/#">Tabela 3</router-link>
                        </div>
                    </li>
                    <li class="drop-down">
                        <p>Narzędzia</p>
                        <div class="drop-down-content">
                            <router-link to="/#">Narzędzie 11</router-link>
                            <router-link to="/#">Narzędzie 2</router-link>
                            <router-link to="/#">Narzędzie 3</router-link>
                        </div>
                    </li>
                </ul>
            </div>
            <div class="login">
                <div class="drop-down">
                    <p>Konto</p>
                    <div class="drop-down-content">
                        <router-link to="/#">Zaloguj</router-link>
                        <router-link to="/#">Zarejestruj</router-link>
                    </div>
                </div>
            </div>
        </nav>
    </div>
</template>

<script>
export default {
    name: 'NavBar',
};
</script>
<style scoped>
nav {
    z-index: 2;
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    line-height: 54px;
    padding: 0 1rem;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.25);
    height: 54px;
}

p {
    margin: 0;
}

.drop-down-content {
    display: none;
    position: absolute;
    background-color: #f9f9f9;
    box-shadow: 0px 8px 16px 0px rgba(0, 0, 0, 0.2);
    border: 1px solid rgba(0, 0, 0, 0.25);
    border-top: 0;
    left: -0.25rem;
    padding: 0 0.25rem;
    white-space: nowrap;
}

.drop-down-content a {
    display: block;
    width: 100%;
    padding: 0 0.25rem;
    text-align: center;
}

.drop-down:hover .drop-down-content {
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
}

.logo a {
    /* position: ; */
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

.logo a::after {
    position: absolute;
    content: '';
    top: 3px;
    left: 1rem;
    width: 100%;
    height: 48px;
    border-right: rgba(0, 0, 0, 0.25) 1px solid;
}

.nav-links ul {
    display: flex;
    justify-content: center;
    align-items: center;
    flex-direction: row;
    flex: 1;
    margin: 0;
    padding: 0;
}

.nav-links ul li {
    list-style: none;
}

.nav-links ul li .nav-link {
    padding: 0 0.75rem;
}

.nav-links a:hover,
.login a:hover {
    background-color: rgba(0, 0, 0, 0.25);
}

.nav-links ul li a,
.drop-down {
    /* padding: 0 0.75rem; */
    margin: 0 1rem;
    display: block;
}

a {
    color: #000;
    text-decoration: none;
}

.login {
    position: relative;
}

@keyframes rotate {
    0% {
        transform: rotate(0deg);
    }

    100% {
        transform: rotate(180deg);
    }
}

@keyframes rotate2 {
    100% {
        transform: rotate(0deg);
    }

    0% {
        transform: rotate(180deg);
    }
}

.login::before {
    position: absolute;
    content: '';
    left: -1rem;
    top: 3px;
    height: 48px;
    background-image: url('@/assets/user.svg');
    background-size: contain;
    background-repeat: no-repeat;
    background-position: center;
    aspect-ratio: 1/1;
    background-size: 1.5rem;
    border-left: rgba(0, 0, 0, 0.25) 1px solid;

}

.drop-down {
    padding: 0 0.75rem;
    cursor: pointer;
    position: relative;
    display: inline-block;
}

.drop-down::after {
    content: '';
    position: absolute;
    background-image: url('@/assets/arrow.svg');
    background-size: contain;
    background-repeat: no-repeat;
    background-position: center;
    aspect-ratio: 1/1;
    top: 40%;
    right: -0.5rem;
    transform: translateY(-50%);
    height: 0.75rem;
    width: 0.75rem;
    animation: rotate2 0.25s forwards ease-in-out;

}

.drop-down:hover::after {
    animation: rotate 0.25s forwards ease-in-out;
}
</style>