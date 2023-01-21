<template>
    <div class="rna-input">
        <div class="rna-input__header">
            <img src="@/assets/settings.svg" alt="settings" class="rna-input__settings-icon" role="button" />
            <div class="rna-input__header__text">
                <span>
                    <h2>Analiza RNA/DNA</h2>
                    <h4>Narzędzie służące do dogłębnej analizy łańcuchów DNA/RNA </h4>
                </span>
            </div>
            <div class="rna-input__header__gradient"></div>
        </div>
        <div class="rna-input__form">
            <div class="rna-input__form__data-type">
                <p>Typ danych</p>
                <select name="type" id="type" v-model="fileType">
                    <option value="text">Tekst</option>
                    <option value="file">Plik</option>
                </select>
            </div>
            <form action="#" method="POST">
                <textarea placeholder="Wpisz RNA lub DNA..." v-if="fileType == 'text'"></textarea>
                <div v-if="fileType == 'file'">
                    <input type="file" />
                </div>
            </form>
            <div class="rna-input__buttons">
                <div class="rna-input__button">Analizuj</div>
                <div class="rna-input__button">Wyczyść</div>
            </div>
        </div>
    </div>
</template>

<script>
export default {
    name: 'RnaInput',
    data() {
        return {
            rnaSequence: '',
            fileType: 'text',
        }
    },
    methods: {
        submit() {
            this.$emit('submit', this.rnaSequence)
        }
    }
}
</script>

<style scoped>
form div {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: center;
    border: 1px solid rgba(0, 0, 0, 0.5);
    height: 206px;
}

.rna-input__form form {
    position: relative;
}

.rna-input__buttons {
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: center;
}

.rna-input__button::after {
    content: '';
    position: absolute;
    top: 50%;
    transform: translate(0, -50%);
    right: 0.75rem;
    width: 1.25rem;
    height: 1.25rem;
    background-image: url('@/assets/accept.svg');
    background-size: contain;
    background-repeat: no-repeat;
}

.rna-input__button:first-child:hover {
    animation: pulse 1s ease-in-out infinite alternate;
}

.rna-input__button:last-child:hover {
    animation: pulse2 1s ease-in-out infinite alternate;
}


@keyframes pulse {
    from {
        box-shadow: 0 0 5px 5px rgba(0, 255, 0, 0.15);
    }

    to {
        box-shadow: 0 0 5px 5px rgba(0, 255, 0, 0.45);
    }
}

@keyframes pulse2 {
    from {
        box-shadow: 0 0 5px 5px rgba(255, 0, 0, 0.15);
    }

    to {
        box-shadow: 0 0 5px 5px rgba(255, 0, 0, 0.45);
    }
}


.rna-input__button:last-child:hover {
    box-shadow: 0 0 5px 5px, rgba(255, 0, 0, 0.1);
}

.rna-input__button:last-child::after {
    background-image: url('@/assets/cancel.svg');
}

.rna-input__button {
    position: relative;
    background-color: #fff;
    border: 1.5px solid rgba(0, 0, 0, 0.5);
    padding: 0.5rem 2.5rem 0.5rem 1.5rem;
    margin: 1rem;
    cursor: pointer;
}

.rna-input__form__data-type {
    padding: 0.5rem;
    display: flex;
    flex-direction: row;
    align-items: center;
    justify-content: flex-end;
}

.rna-input__form__data-type select {
    margin-left: 0.5rem;
    height: fit-content;
    background-color: #fff;
    border: 1.5px solid rgba(0, 0, 0, 0.5);
    appearance: none;
    background-image: url('@/assets/arrow.svg');
    background-size: 0.5rem;
    background-repeat: no-repeat;
    background-position-x: 95%;
    background-position-y: 50%;
    padding: 0.4rem;
    padding-right: 2.4rem;
    cursor: pointer;
}

option:hover {
    cursor: pointer;
}

.rna-input__form__data-type select:hover {
    filter: brightness(0.9);
}

.rna-input__settings-icon {
    width: 1.5rem;
    height: 1.5rem;
    position: absolute;
    top: 0.5rem;
    right: 0.5rem;
    cursor: pointer;
}

.rna-input__settings-icon:hover {
    transition: 0.5s;
    filter: brightness(0);
}

.rna-input {
    border: 1.5px solid rgba(0, 0, 0, 0.5);
    width: 70%;
}

h2,
h4 {
    margin: 0;
}

h2 {
    font-weight: 400;
}

h4 {
    font-weight: 200;
}

.rna-input__header {
    position: relative;
    display: flex;
    flex-direction: row;
}

.rna-input__header__text {
    padding: 0.5rem;
    background-color: #5F5CFF;
    color: white;
}

.rna-input__header__gradient {
    background: rgb(95, 92, 255);
    background: linear-gradient(90deg, rgba(95, 92, 255, 1) 0%, rgba(255, 255, 255, 1) 100%);
    flex: 1;
}

textarea {
    overflow-y: auto;
    width: calc(100% - 0.5rem);
    resize: none;
    height: 200px;
    box-shadow: inset 0px 0px 10px 0px rgba(0, 0, 0, 0.5);
    margin: 0;
}

.rna-input__form {
    padding: 0 1.5rem 0rem 1.5rem;
}
</style>