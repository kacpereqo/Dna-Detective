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

            axios.post('http://localhost:8000/login', {
                username: email,
                password: password
            }, {
                headers: {
                    "Content-Type": "application/x-www-form-urlencoded"
                }
            }).then(res => {
                const jwt = res.data.access_token;
                localStorage.setItem('jwt', jwt);
                this.$store.commit('setUser', jwt);
                this.$router.push('/');
            }).catch(err => {
                this.error = err.response.data.detail;
            })
        }
    }
}
</script>

<style scoped>
.login-wrapper {
    border: 1px solid var(--accent-color-dark);
    border-radius: 5px;
    width: 300px;
    height: 350px;
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
    text-align: left;
    margin-bottom: 20px;
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