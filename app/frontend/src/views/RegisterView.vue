<template>
    <div class="wrapper">
        <div class="login-wrapper">
            <div class="header">
                <h1>Rejestracja</h1>
            </div>
            <form>
                <input type="text" placeholder="E-mail" />
                <input type="password" placeholder="Hasło" />
                <button type="submit">Zaloguj</button>
            </form>
            <div class="info">
                <p>Masz juz konto? <router-link to="/login">Zaloguj się</router-link></p>
            </div>
        </div>
    </div>
</template>

<script>
export default {
    name: 'Login',
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
    /* min-width: 20vw; */
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
    border: 1px solid var(--accent-color);
    border-radius: 5px;
    background-color: transparent;
}

.login-wrapper button {
    padding: 10px;
    margin: 10px 0;
    font-size: 1rem;
    border: 1px solid var(--accent-color);
    border-radius: 5px;
    background-color: transparent;
    color: var(--text-color);
    cursor: pointer;
}

.info {
    font-size: 0.9rem;
}
</style>