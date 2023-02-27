<template>
    <div class="wrapper">
        <div class="profile-wrapper">
            <div class="sequences">
                <ul>
                    <li>
                        <h2>Historia wykonanych translacji</h2>
                    </li>
                    <li @click="showPreview" v-for="(value, key) in data.sequences" :key="key" class="sequence">
                        <div><span>{{ value }}</span>
                            <router-link :to="`/translations/${value}`">
                                <img src="@/assets/link.svg">
                            </router-link>
                        </div>
                    </li>
                </ul>
            </div>
            <div class="preview">
                <h2>Podgląd</h2>
                <div class="preview-content">
                    <span v-if="selected.length == 0">Klilnij na sekwencję aby zobaczyć jej podgląd</span>
                    <span v-else>{{ selected }}</span>
                </div>
            </div>
        </div>
    </div>
</template>

<script>
import axios from 'axios';
import { useMeta } from 'vue-meta'

export default {
    name: "ProfileView",
    setup() {
        useMeta({
            title: 'Profile',

        })
    },
    data() {
        return {
            data: {},
            selected: "",
        }
    },
    mounted() {
        axios.get('http://127.0.0.1:8000/api/profile').then(res => {
            this.data = res.data
        }).catch(err => {
            console.log(err)
        })
    },
    methods: {
        showPreview() {
            axios.get(`http://127.0.0.1:8000/api/63fb502c37de766e5bd6d93b/translate?is_reversed=false&is_forward=true`).then(res => {
                this.selected = res.data
            }).catch(err => {
                console.log(err)
            })
        }
    }
}

</script>

<style scoped>
.wrapper {
    align-items: stretch;
    align-content: stretch;
}

.profile-wrapper {
    flex: 1;
    display: flex;
    padding: 2rem;
}

.sequences {
    overflow: auto;
    border-radius: 0.5rem;
    border: 1px solid var(--accent-color-dark);
    flex: 2;
    margin-right: 3rem;
}

.preview {
    display: flex;
    flex-direction: column;
    align-items: stretch;
    overflow: hidden;
    border-radius: 0.5rem;
    border: 1px solid var(--accent-color-dark);
    flex: 5;
}

.preview-content {
    flex: 1;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
}

.sequence {
    word-break: break-all;
    padding: 1rem;
    border-bottom: 1px solid var(--accent-color-dark);
}

.sequences ul {
    list-style: none;
    padding: 0;
    margin: 0;
}


h2 {
    color: white;
    padding: 0.5rem;
    margin: 0;
    font-size: 1.3rem;
    background-color: var(--main-color);
}

.sequences li:hover {
    background-color: var(--accent-color-light);
    cursor: pointer;
}

li div {
    display: flex;
    justify-content: space-between;
    align-items: center;
}

img {
    width: 1.5rem;
    height: 1.5rem;
    filter: var(--icon-filter);
}

img:hover {
    cursor: pointer;
    box-shadow: 0 0 0px 4px var(--accent-color-dark);
    border-radius: 100%;
}
</style>