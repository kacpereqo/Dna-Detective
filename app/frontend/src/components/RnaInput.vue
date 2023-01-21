<template>
    <div class="rna-input">
        <div class="rna-input__header">
            <img <div class="rna-input__header__text">
            <span>
                <h2>Analiza RNA/DNA</h2>
                <h4>Narzędzie służące do dogłębnej analizy łańcuchów DNA/RNA </h4>
            </span>
        </div>
        <div class="rna-input__header__gradient"></div>
    </div>
    <div class="rna-input__form">
        <form action="#" method="POST">
            <textarea placeholder="Wpisz RNA lub DNA">

            </textarea>
        </form>
    </div>
    </div>
</template>

<script>
export default {
    name: 'RnaInput',
    data() {
        return {
            rnaSequence: ''
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
}

.rna-input__form {
    padding: 1.5rem;
}
</style>