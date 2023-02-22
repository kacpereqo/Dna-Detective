<template>
    <div class="sidebar-wrapper" :style="{ 'left': showMenu ? '0' : '-217px' }">
        <div class="show-menu" @click="showMenu = !showMenu"><span
                :style="{ '--display': showMenu ? 'none' : 'block' }"></span></div>
        <div class="items" :style="{ 'position': position }">
            <SequenceName />
            <input type="text" placeholder="Szukaj" id="search" v-model="search" />
            <ul>
                <li v-for="item in list" @click="changeComponent(item.value)" class="parent-li">
                    <p class="heading">{{
                        item.text
                    }}</p>
                    <ul>
                        <li class="sub-li" v-for="item, index in item.nested">
                            <p @click="scrollToContent(item.value)" :style="{ '--i': index }">{{ item.text }}</p>
                        </li>

                    </ul>
                </li>
            </ul>
        </div>

    </div>
</template>

<script>
import SequenceName from '@/components/navigation/SequenceName.vue';

export default {
    name: 'Sidebar',

    components: {
        SequenceName,
    },

    data() {
        return {
            showMenu: false,
            position: 'inherit',

            list: [
                {
                    text: 'Wizualizacja', value: "visualization", nested: [
                        { text: 'Obraz', value: "" }
                    ]
                },
                {
                    text: 'Ładunek', value: "charge", nested: [
                        { text: 'Ładunek od ph', value: "" },
                        { text: 'Punkt izoelektryczny', value: "" }
                    ]
                },
                {
                    text: 'Hydrofobowość', value: "hydro", nested: [
                        { text: 'Wykres Hydrofobowości', value: "" },
                        { text: 'Średnia Hydrofobowość', value: "" }
                    ]
                },
                {
                    text: 'Właściwości', value: "propeties", nested: [
                        { text: 'Masa', value: "" }
                    ]
                },
            ],
            search: '',
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
                this.position = 'fixed';
            } else {
                this.position = 'inherit';
            }
        },
        changeComponent(componentName) {
            this.$emit('changeComponent', componentName);
        },
        scrollToContent(element) {
            this.$emit('scrollToContent', element);
        }
    }

}
</script>

<style scoped>
#search {
    background: transparent;
    width: calc(100% - 1.5rem);
    padding: 0.5rem;
    border: var(--accent-color) 1px solid;
    border-radius: 0.25rem;
    margin-bottom: 1.5rem;
}

.sidebar-wrapper {
    overflow: hidden;
    background-color: var(--background-color);
    margin-top: 2px;
    flex-shrink: 0;
    width: 13.5rem;
    min-height: calc(100vh - 64px);
    border-right: var(--accent-color) 1px solid;
}

.items {
    left: 0;
    top: 0;
    width: 12rem;
    /* width: inherit; */
    display: flex;
    flex-direction: column;
    align-items: center;
    overflow: scroll;
    overflow-x: hidden;
    scrollbar-width: thin;
    justify-content: flex-start;
    height: calc(100vh - 56px);
    padding: 0.75rem;
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
    margin: 0.25rem 0;
}


ul li {
    margin: 0;
    padding-left: 0.25rem;
    border-left: 1px solid var(--accent-color);
}

ul li ul li {
    border-left: none !important;
}


li p {
    display: block;
    white-space: nowrap;
    font-size: 0.9rem;
    display: block;
    padding: 0.5rem 0.25rem;
    margin: 0;
}

li p:hover {
    cursor: pointer;
}

p:hover {
    text-decoration: underline;
}

.heading:hover::after {
    content: "";
    position: absolute;
    top: 0;
    left: -2px;
    height: 36px;
    border-right: 2px solid;
}

.sub-li:hover p::after {
    content: "";
    position: absolute;
    top: calc(var(--i) * 36px + 41px);
    left: -2px;
    height: 36px;
    border-right: 2px solid;
}

.parent-li li:hover .heading::after {
    display: none;
}

li {
    margin: 0 0.5rem;
}

.heading {
    font-size: 1.15rem;
}

.parent-li {
    position: relative;
    margin-bottom: 0.5rem;
}

@media (max-width: 960px) {
    .sidebar-wrapper {
        position: absolute;
        z-index: 1;
        left: -217px;
    }

    .show-menu {
        display: block !important;
    }
}

.show-menu {
    display: none;
    position: absolute;
    ;
    right: calc(-2.5rem - 2px);
    top: 0.25rem;
    background-color: var(--background-color);
    border: 1px solid var(--accent-color);
    border-left: none;
    border-radius: 0 0.25rem 0.25rem 0;
    width: 2.5rem;
    height: 2.5rem;
}

.show-menu span {
    display: var(--display);
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    border-top: 2px solid var(--background-color);
    width: 2rem;
}

.show-menu:hover {
    cursor: pointer;
}
</style>