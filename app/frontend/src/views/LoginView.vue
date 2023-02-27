<template>
    <div class="wrapper">
        <div class="login-wrapper">
            <div class="header">
                <h1>Logowanie</h1>
            </div>
            <form>
                <input type="text" ref="email" placeholder="E-mail" />
                <input type="password" ref="password" placeholder="Hasło" />
                <div class="error-message">{{ error }}</div>
                <button type="submit" @click="submit">Zaloguj</button>
            </form>
            <div class="info">
                <p>Nie masz konta? <router-link to="/register">Zarejestruj się</router-link></p>
            </div>
        </div>
    </div>
</template>

<script>
import axios from 'axios';
import { useMeta } from 'vue-meta'

export default {
    name: 'Login',
    setup() {
        useMeta({
            title: 'Logowanie',
        })
    },
    data() {
        return {
            error: ''
        }
    },

    methods: {
        submit() {
            this.error = '';
            const email = this.$refs.email.value;
            const password = this.$refs.password.value;

            axios.post('http://localhost:8000/token', {
                username: email,
                password: password
            }, {
                headers: {
                    "Content-Type": "application/x-www-form-urlencoded"
                }
            }).then(res => {
                this.$store.commit("setEmail", res.data.email)
                console.log(res.data.preferences)
                this.$store.commit('toogleTheme', res.data.preferences.theme)
                this.$i18n.locale = res.data.preferences.lang;
                localStorage.setItem('locale', res.data.preferences.lang);
                this.$store.commit('setFontSize', res.data.preferences.fontSize);

                const jwt = res.data.access_token;
                localStorage.setItem('jwt', jwt);
                this.$store.commit('setUser', jwt);
                this.$router.push('/');
            }).catch(err => {
                console.log(err)
            })
        }
    }
}
</script>

<style scoped>
a {
    text-decoration: underline;
}

.login-wrapper {
    overflow: hidden;
    border: 1px solid var(--accent-color-dark);
    border-radius: 5px;
    width: 20rem;
    height: 25rem;
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    align-items: center;
}

.login-wrapper form {
    width: 100%;
    padding: 0 10px;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
}

.header {
    overflow: hidden;
    text-align: left;
    margin-bottom: 20px;
    height: 2.5rem;
    background-color: var(--main-color);
    color: white;
    width: 100%;
}

.header h1 {
    font-size: 1.5rem;
    font-weight: normal;
    margin: 0;
    padding: 10px;
}

.login-wrapper input {
    width: 70%;
    padding: 10px;
    font-size: 1rem;
    margin: 10px 0;
    background-color: transparent;
    border: 1px solid var(--accent-color);
    border-radius: 5px;
}

.login-wrapper button {
    padding: 10px;
    margin: 10px 0;
    font-size: 1rem;
    background-color: transparent;
    border: 1px solid var(--accent-color);
    border-radius: 5px;
    color: var(--text-color);
    cursor: pointer;
    transition: 0.3s;
}

.login-wrapper button:hover {
    background-color: var(--accent-color-light);
}

.info {
    font-size: 0.9rem;
}

.error-message {
    color: red;
}
</style>