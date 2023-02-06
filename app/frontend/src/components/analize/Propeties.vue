<template>
    <div class="propeties-wrapper">
        <ChartWrapper :data="weight" :labels="labels" :xTitle="'Aminowkas'" :yTitle="'Masa [Dal]'" />
    </div>
</template>

<script>
import axios from 'axios'
import ChartWrapper from '../ChartWrapper.vue';

export default {
    name: "Propeties",
    data() {
        return {
            weight: [],
            labels: [],
            window: 3,
        };
    },
    methods: {
        getWeight() {
            axios.get(`http://127.0.0.1:8000/api/weight/${this.id}?window=${this.window}`)
                .then(res => {
                    this.weight = res.data.weight;
                    for (let i = this.window; i < this.window + this.weight.length; i++) {
                        this.labels.push(i);
                    }
                });
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getWeight();
    },
    components: { ChartWrapper }
}

</script>

<style scoped>
.propeties-wrapper {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
}
</style>