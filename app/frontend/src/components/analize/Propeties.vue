<template>
    <div class="propeties-wrapper">
        <ChartWrapper v-if="loaded" :xData="xData" :yData="yData" :labels="labels" :wholeNumbers="true" />
    </div>
</template>

<script>
import axios from 'axios'
import ChartWrapper from '../chart/ChartWrapper.vue';

export default {
    name: "Propeties",
    data() {
        return {
            yData: [],
            xData: [],
            window: 3,
            loaded: false,
        };
    },
    methods: {
        getWeight() {
            axios.get(`http://127.0.0.1:8000/api/weight/${this.id}?window=${this.window}`)
                .then(res => {
                    this.yData = res.data.weight;
                    for (let i = this.window; i < this.window + this.yData.length; i++) {
                        this.xData.push(i);
                    }
                    this.loaded = true;
                });
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getWeight();
    },
    mounted() {
        this.labels = this.$t('charts.mass');
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