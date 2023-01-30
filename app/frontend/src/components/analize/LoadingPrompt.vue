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
import LoadingIcon from '../LoadingIcon.vue';

export default {
    name: 'LoadingPrompt',
    components: {
        LoadingIcon,
    },
    data() {
        return {
            dots: 0,
            fact: "W 2001 roku, w wyniku badań genetycznych, odkryto, że ludzie mają 69 chromosomów, a nie 48 jak się wcześniej myślało."
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
        }
    },
    methods: {
        addDots() {
            this.dotsInterval = setInterval(() => {
                if (this.dots < 3) {
                    this.dots += 1;
                } else {
                    this.dots = 0;
                }
            }, 500);
        },
        getFact() {
            this.factInterval = setInterval(() => {
                axios.get('http://127.0.0.1:8000/fact')
                    .then(response => {
                        this.fact = response.data.fact;
                    })
                    .catch(error => {
                        console.log(error);
                    });
            }, 10000);

        }
    },
    created() {
        this.addDots();
        if (this.display_facts) {
            this.getFact();
        }
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
    font-size: 0.8rem;
    text-align: center;
    width: 40%;
}
</style>