<template>
    <div class="sidebar-wrapper">
        <div class="sticky" :style="{ 'position': isFixed }">
            <div class="items">
                <input type="text" placeholder="Szukaj" id="search" />
                <ul>
                    <li>
                        Wizualizacja
                    </li>
                    <li>
                        Hydrofobowość
                    </li>
                    <li>
                        Ładunek
                    </li>
                </ul>
            </div>
        </div>
    </div>
</template>

<script>

export default {
    name: 'Sidebar',

    data() {
        return {
            isFixed: 'inherit',
        };
    },

    mounted() {
        window.addEventListener('scroll', this.handleScroll);
    },

    beforeDestroy() {
        window.removeEventListener('scroll', this.handleScroll);
    },

    methods: {
        handleScroll() {
            const sidebar = document.querySelector('.sidebar-wrapper');
            const sidebarTop = sidebar.getBoundingClientRect().top;

            if (sidebarTop <= 0) {
                this.isFixed = 'fixed';
            } else {
                this.isFixed = 'inherit';
            }
        },
    },

}
</script>

<style scoped>
#search {
    width: calc(100% - 1rem);
    padding: 0.5rem;
    border: rgba(0, 0, 0, 0.25) 1px solid;
    border-radius: 0.25rem;
    margin-bottom: 0.5rem;
}

.sidebar-wrapper {
    width: 216px;
    border-right: rgba(0, 0, 0, 0.25) 1px solid;
}

.items {
    display: flex;
    flex-direction: column;
    align-items: flex-start;
    justify-content: flex-start;
    height: 100%;
    padding: 1rem;
}

.sticky {
    width: inherit;
    position: sticky;
    top: 0;
    overflow: auto;
    align-self: flex-start;
}

ul {
    list-style-type: none;
    padding: 0;
    margin: 0 0 0 0.5rem;
}

li {
    margin: 0.5rem 0;
}
</style>