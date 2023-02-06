<template>
    <div class="sidebar-wrapper" :class="{ 'fixed': position }">
        <div class="items">
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
import SequenceName from '@/components/SequenceName.vue';

export default {
    name: 'Sidebar',

    components: {
        SequenceName,
    },

    data() {
        return {
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
                this.isFixed = 'fixed';
            } else {
                this.isFixed = 'inherit';
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
    margin-top: 2px;
    flex-shrink: 0;
    width: 216px;
    border-right: var(--accent-color) 1px solid;
}

.items {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: flex-start;
    height: 100%;
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
</style>