<template>
    <div class="rna-input">
        <div class="rna-input__header">
            <img src="@/assets/settings.svg" alt="settings" class="rna-input__settings-icon" role="button" />
            <div class="rna-input__header__text">
                <span>
                    <h1>Analiza RNA/DNA</h1>
                    <h2>Narzędzie służące do dogłębnej analizy łańcuchów DNA/RNA </h2>
                </span>
            </div>
        </div>
        <div class="rna-input__form">
            <div class="rna-input__form__data-type">
                <div class="messages">
                    <p v-if="errorMessage">{{ errorMessage }} <img src="@/assets/error.svg"></p>
                </div>
                <div class="data-type">
                    <p>Typ danych</p>
                    <select name="type" id="type" v-model="fileType">
                        <option value="text">Tekst</option>
                        <option value="file">Plik</option>
                    </select>
                </div>
            </div>
            <form action="#" method="POST">
                <div v-if="fileType == 'text'">
                    <textarea placeholder="Wpisz RNA lub DNA..." v-model="rnaSequence" spellcheck="false"></textarea>
                </div>
                <div v-if="fileType == 'file'">
                    <DropFile />
                </div>
            </form>
            <div class="rna-input__buttons">
                <div class="rna-input__button" role="button" @click="submit"
                    style="--first-color:rgba(0,255,0,0.25); --second-color:rgba(0,255,0,0.4);">Analizuj</div>
                <div class="rna-input__button" role="button" @click="clear"
                    style="--first-color:rgba(255,0,0,0.25); --second-color:rgba(255,0,0,0.4);">Wyczyść</div>
            </div>
        </div>
    </div>
</template>

<script>
import DropFile from '@/components/rnainput/DropFile.vue';
import axios from 'axios'

export default {
    name: 'RnaInput',
    components: {
        DropFile
    },
    data() {
        return {
            rnaSequence: '',
            fileType: 'text',
            errorMessage: '',
        }
    },
    methods: {
        submit() {
            this.validate().then(() => {
                axios.post('http://127.0.0.1:8000/api/sequence', {
                    sequence: this.rnaSequence
                }).then(res => {
                    this.$router.push({
                        name: 'translations',
                        params: {
                            id: res.data.id
                        }
                    })

                }).catch(err => {
                    if (err.response) {
                        this.errorMessage = err.response.data.detail
                    }
                })
            }).catch(err => {
                this.errorMessage = err.message
            })
        },
        clear() {
            this.rnaSequence = ''
        },
        validate() {
            return new Promise((resolve, reject) => {
                if (this.rnaSequence.length < 1) {
                    reject(new Error('Pole nie może być puste'))
                }

                if (this.fileType == 'text') {
                    if (this.rnaSequence.length < 3) {
                        reject(new Error('Sekwencja musi zawierać przynajmniej 3 znaki'))
                    }

                    if (!this.rnaSequence.match(/^[AUGCT]+$/i)) {
                        reject(new Error('Sekwencja może zawierać tylko znaki A, U, G, C, T'))
                    }
                }

                resolve()
            })
        }

    }
}
</script>

<style scoped lang="scss">
form div {
    border: 1px solid var(--accent-color-dark);
    height: 200px;
}

@media screen and (max-width: 1024px) {
    .rna-input {
        width: 100% !important;
    }
}

.rna-input__buttons {
    display: flex;
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
    filter: var(--icon-filter);
    background-size: contain;
}

.rna-input__button:hover {
    background-color: var(--accent-color-light);
    animation: pulse ease-in-out .7s infinite alternate;
}

@keyframes pulse {
    from {
        box-shadow: 0 0 0 0 var(--first-color);
    }

    to {
        box-shadow: 0 0 10px 3px var(--second-color);
    }
}

.rna-input__button:last-child::after {
    background-image: url('@/assets/cancel.svg');
}

.rna-input__button {
    position: relative;
    border: 1.5px solid var(--accent-color-dark);
    padding: 0.5rem 2.5rem 0.5rem 1.5rem;
    margin: 1rem;
    cursor: pointer;
}

.rna-input__form__data-type {
    margin: 0.5rem 0;
    display: flex;
    align-items: center;
    justify-content: center;
}

.data-type {
    display: flex;
    align-items: center;
}

.messages {
    flex: 3;
    color: rgb(255, 0, 0);

}

.messages p {
    width: fit-content;
    position: relative;
}

.messages img {
    bottom: 0;
    right: -2rem;
    position: absolute;
    filter: brightness(0) saturate(100%) invert(16%) sepia(75%) saturate(6903%) hue-rotate(357deg) brightness(107%) contrast(115%);
}

.rna-input__form__data-type select {
    margin-left: 0.5rem;
    height: fit-content;
    background-color: transparent;
    border: 1.5px solid var(--accent-color-dark);
    appearance: none;
    background-image: url('@/assets/arrow.svg');
    filter: var(--icon-filter);
    background-size: 0.5rem;
    background-repeat: no-repeat;
    background-position-x: 95%;
    background-position-y: 50%;
    padding: 0.4rem;
    padding-right: 3rem;
    cursor: pointer;
}

select option:hover {
    cursor: pointer;
}

.rna-input__form__data-type select:hover {
    filter: brightness(0.85);
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
    border: 1px solid var(--accent-color-dark);
    width: 70%;
}

h1,
h2 {
    margin: 0;
}

h1 {
    font-size: 1.4rem;
    font-weight: 400;
}

h2 {
    font-weight: 200;
    font-size: 1rem;
}

.rna-input__header {
    position: relative;
    display: flex;
    color: white;
}

.rna-input__header__text {
    padding: 0.5rem;
    background-color: var(--main-color);
    width: 100%;
}



textarea {
    border: 0;
    overflow-y: auto;
    width: 100%;
    resize: none;
    padding: 0;
    height: 100%;
    margin: 0;
    border: 1px solid var(--accent-color-dark);
    background-color: transparent;
}

::placeholder {
    color: var(--accent-color-dark);
}

.rna-input__form {
    padding: 0 1.5rem 0rem 1.5rem;
}
</style>