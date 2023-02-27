<template>
    <div class="loading-prompt-wrapper">
        <LoadingIcon />
        <div class="text">{{ text }} <span class="dots">{{ '.'.repeat(dots) }}</span>
        </div>
        <div class="facts" v-if="display_facts">
            {{ fact }}
        </div>
    </div>
</template>

<script>
import axios from 'axios';
import LoadingIcon from './LoadingIcon.vue';

export default {
    name: 'LoadingPrompt',
    components: {
        LoadingIcon,
    },
    data() {
        return {
            dots: 0,
            factInterval: null,
            dotsInterval: null,
            fact: "W 2001 roku, w wyniku badań genetycznych, odkryto, że ludzie mają 48 chromosomów, a nie 69 jak się wcześniej myślano."
        }
    },
    props: {
        text: {
            type: String,
            required: false,
            default: 'Ładowanie'
        },
        display_facts: {
            type: Boolean,
            required: false,
            default: false
        },
        loaded: {
            type: Boolean,
            required: false,
            default: false
        }
    },
    methods: {
        addDots() {
            if (this.dots < 3) {
                this.dots += 1;
            } else {
                this.dots = 0;
            }
        },
        getFact() {
            axios.get('http://127.0.0.1:8000/fact')
                .then(response => {
                    this.fact = response.data.fact;
                })
        },
        clearIntervals() {
            clearInterval(this.dotsInterval);
            clearInterval(this.factInterval);
        }
    },
    created() {
        this.dotsInterval = setInterval(this.addDots, 500);
        this.factInterval = setInterval(this.getFact, 1000);
    },
    beforeUnmount() {
        this.clearIntervals();
    },
}
</script>

<style scoped>
.loading-prompt-wrapper {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    height: 100%;
}

.text {
    font-size: 1.5rem;
    margin-top: 2rem;
    letter-spacing: 0.1rem;
    font-weight: bolder;
}

.dots {
    position: absolute;
}

.facts {
    margin-top: 0.5rem;
    font-size: 0.9rem;
    text-align: center;
    width: 40%;
}
</style>