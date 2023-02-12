<template>
    <hr :style="{ '--word': `'${word}'`, 'border-top': line ? '3px solid var(--accent-color)' : '0' }">
</template>

<script>
export default {
    name: 'Line',
    props: {
        word: {
            type: String,
            default: ''
        },
        line: {
            type: Boolean,
            default: true
        }
    }
}
</script>

<style scoped>
hr {
    border: none;
    color: var(--text-color);
    width: 100%;
    overflow: visible;
    text-align: center;
    height: 5px;
}

hr:after {
    content: var(--word);
    background-color: var(--background-color);
    padding: 0 0.5rem;
    font-size: 1.5rem;
    position: relative;
    top: -1.25rem;
}
</style>