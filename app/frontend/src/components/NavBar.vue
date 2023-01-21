<template>
    <div class="nav-wrapper">
        <nav>
            <div class="logo">
                <router-link to="/">
                </router-link>
            </div>
            <div class="nav-links">
                <ul>
                    <li><router-link to="/">Strona główna</router-link></li>
                    <li><router-link to="/#">Analiza</router-link></li>
                    <li><router-link to="/#" class="drop-down">Tablice</router-link></li>
                    <li><router-link to="/#" class="drop-down">Narzędzia</router-link></li>
                </ul>
            </div>
            <div class="login">
                <a href="#">Login</a>
            </div>
        </nav>
    </div>
</template>

<script>
</script>

<style scoped>
nav {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    line-height: 54px;
    padding: 0 1rem;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.25);
    height: 54px;
}

.logo a {
    position: absolute;
    display: block;
    height: 48px;
    background-image: url('../assets/logo.png');
    background-size: contain;
    background-repeat: no-repeat;
    background-position: center;
    aspect-ratio: 4/1;
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

@media screen and (max-width: 768px) {}

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
    margin: 0 1rem;
}

a {
    color: #000;
    text-decoration: none;
}

.login {
    position: relative;
}

.login::before {
    position: absolute;
    content: '';
    left: -1rem;
    top: 3px;
    height: 48px;
    border-left: rgba(0, 0, 0, 0.25) 1px solid;

}

.drop-down {
    position: relative;
}

.drop-down::after {
    content: '';
    position: absolute;
    background-image: url('@/assets/arrow.svg');
    background-size: contain;
    background-repeat: no-repeat;
    background-position: center;
    aspect-ratio: 1/1;
    top: 50%;
    right: -1rem;
    transform: translateY(-50%);
    height: 0.75rem;
    width: 0.75rem;

}
</style>