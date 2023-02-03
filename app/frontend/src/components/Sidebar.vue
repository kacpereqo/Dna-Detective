<template>
    <div class="sidebar-wrapper" :class="{ 'fixed': position }">
        <div class="items">
            <input type="text" placeholder="Szukaj" id="search" v-model="search" />
            <ul>
                <li v-for="item in list" @click="changeComponent(item.value)"><span class="heading">{{
                    item.text
                }}</span>
                    <ul v-for="item in item.nested">
                        <li><span @click="scrollToContent(item.value)">{{ item.text }}</span></li>
                    </ul>
                </li>
            </ul>
        </div>
    </div>
</template>

<script>

export default {
    name: 'Sidebar',

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
    width: calc(100% - 1rem);
    padding: 0.5rem;
    border: rgba(0, 0, 0, 0.25) 1px solid;
    border-radius: 0.25rem;
    margin-bottom: 0.5rem;
}

.sidebar-wrapper {
    flex-shrink: 0;
    width: 216px;
    border-right: rgba(0, 0, 0, 0.25) 1px solid;
}

.items {
    display: flex;
    flex-direction: column;
    align-items: flex-start;
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
    margin: 0.65rem 0;
}


ul li {
    padding-left: 0.25rem;
    border-left: 1px solid rgba(0, 0, 0, 0.25);
}

ul li ul li {
    border-left: none !important;
}


li span {
    white-space: nowrap;
    font-size: 0.9rem;
    display: block;
    padding: 0.15rem 0;
    border-radius: 0.1rem;
}

span:hover {
    cursor: pointer;
    font-weight: bold;
}

.heading {
    font-size: 1.15rem;
}
</style>