<template>
    <div class="wrapper">
        <div class="login-wrapper">
            <div class="header">
                <h1>Rejestracja</h1>
            </div>
            <form @submit="submit">
                <input type="text" ref="login" placeholder="Login" />
                <input type="text" placeholder="E-mail" ref="email" />
                <input type="password" placeholder="Hasło" ref="password" />
                <div class="error">{{ error }}</div>
                <button type="submit">Zarejestruj</button>
            </form>
            <div class="info">
                <p>Masz juz konto? <router-link to="/login">Zaloguj się</router-link></p>
            </div>
        </div>
    </div>
</template>

<script>
import axios from 'axios';
import { useMeta } from 'vue-meta'

export default {
    name: 'Register',
    setup() {
        useMeta({
            title: 'Rejestracja',
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

            axios.post('http://localhost:8000/register', {
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
    padding: 0.5rem;
    height: 2.5rem;
    overflow: hidden;
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
    height: 2rem;
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