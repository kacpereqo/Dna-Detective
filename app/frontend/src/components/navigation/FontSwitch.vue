<template>
    <div class="switch-wrapper">
        <img src="@/assets/text_decrease.svg" @click="changeFont(-2)">
        <img src="@/assets/text_increase.svg" @click="changeFont(2)">
    </div>
</template>
<script>
import axios from 'axios';
export default {
    name: "FontSwitch",

    methods: {
        changeFont(size) {
            const root = document.documentElement;
            const currentSize = parseInt(window.getComputedStyle(root).getPropertyValue('font-size').replace('px', ''));
            if (currentSize + size >= 10 && currentSize + size <= 20) {
                root.style.fontSize = currentSize + size + 'px';
                this.updateUserPreference(currentSize + size);
                localStorage.setItem('font-size', currentSize + size + 'px');
            }
        }, updateUserPreference(fontSize) {
            if (this.$store.state.isLogged) {
                axios.post(`http://127.0.0.1:8000/api/preferences?key=fontSize&value=${fontSize}`, {
                })
            }
        }
    }
}
</script>

<style scoped>
.switch-wrapper {
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
}

.switch-wrapper img:first-child {
    margin-right: 0.4rem;
}

.switch-wrapper img {
    width: 1.9rem;
    height: 1.9rem;
    filter: var(--icon-filter);
}
</style>