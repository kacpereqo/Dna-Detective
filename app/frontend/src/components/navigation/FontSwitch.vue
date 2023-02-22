<template>
    <div class="switch-wrapper">
        <img src="@/assets/text_decrease.svg" @click="changeFont(-2)">
        <img src="@/assets/text_increase.svg" @click="changeFont(2)">
    </div>
</template>
<script>
export default {
    name: "FontSwitch",

    methods: {
        changeFont(size) {
            const root = document.documentElement;
            const currentSize = parseInt(window.getComputedStyle(root).getPropertyValue('font-size').replace('px', ''));
            if (currentSize + size >= 10 && currentSize + size <= 20) {
                root.style.fontSize = currentSize + size + 'px';
                localStorage.setItem('font-size', currentSize + size + 'px');
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
    margin: 0 1rem;
    cursor: pointer;
}

.switch-wrapper img {
    width: 2.5rem;
    height: 2.5rem;
    filter: var(--icon-filter);
}
</style>