<template>
    <div class="sequence-name">
        <div class="input-wrapper" @dblclick="editName">
            <input v-model="name" ref="name">
            <img src="../../assets/edit.svg" @click="focusInput">
        </div>
    </div>
</template>

<script>

export default {
    name: 'SequenceName',
    data() {
        return {
            name: 'seq-751931',
        }
    },
    methods: {
        focusInput() {
            this.$refs.name.focus();
        }
    }
}

</script>
aa
<style scoped>
input {
    background-color: transparent;
    text-align: center;
    box-sizing: border-box;
    font-size: 1rem;
    margin: 0;
    font-weight: bold;
    font-family: 'Nunito', sans-serif;
    border: none;
    width: 10rem;
    padding: 0;
}

input:focus {
    outline: none;
    box-shadow: 0px 1px 0px var(--accent-color), 0px 2px 0px var(--accent-color-light);
}

.sequence-name {
    width: 100%;
    margin-bottom: 1rem;
}

.input-wrapper {
    width: inherit;
    position: relative;
    display: flex;
    align-items: center;
    justify-content: center;
}

h1 {
    margin: 0;
    font-size: 1rem;
}

img {
    width: 1.4rem;
    height: 1.4rem;
    cursor: pointer;
    margin-left: 0.25rem;
    filter: var(--icon-filter)
}
</style>