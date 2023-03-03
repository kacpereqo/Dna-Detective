<template>
    <div class="propeties-wrapper">
        <ChartWrapper v-if="loaded.weight" :xData="xData.weight" :yData="yData.weight" :labels="'charts.weight'"
            :wholeNumbers="true" />
        <ChartWrapper v-if="loaded.bulkiness" :xData="xData.bulkiness" :yData="yData.bulkiness" :labels="'charts.bulkiness'"
            :wholeNumbers="true" />
        <ChartWrapper v-if="loaded.recognition" :xData="xData.recognition" :yData="yData.recognition"
            :labels="'charts.recognition'" :wholeNumbers="true" />
    </div>
</template>

<script>
import axios from 'axios'
import ChartWrapper from '../chart/ChartWrapper.vue';

export default {
    name: "Propeties",
    data() {
        return {
            yData: {},
            xData: {},
            window: 3,
            loaded: {
                weight: false,
                bulkiness: false,
                recognition: false,
            },
            labels: {},
        };
    },
    methods: {
        getWeight() {
            axios.post(`http://127.0.0.1:8000/api/weight?window=${this.window}`, {
                frame: this.$store.state.frame,
            })
                .then(res => {
                    this.yData.weight = [];
                    this.xData.weight = [];

                    this.yData.weight = res.data.weight;
                    for (let i = 2; i < 2 + this.yData.weight.length; i++) {
                        this.xData.weight.push(i);
                    }
                    this.loaded.weight = true;
                });
        },
        getBulkiness() {
            axios.post(`http://127.0.0.1:8000/api/bulkiness`, {
                frame: this.$store.state.frame,
            })
                .then(res => {
                    this.yData.bulkiness = [];
                    this.xData.bulkiness = [];

                    for (let x in res.data.bulkiness) {
                        this.yData.bulkiness.push(res.data.bulkiness[x]);
                        this.xData.bulkiness.push(parseInt(x) + 1);
                    }
                    this.loaded.bulkiness = true;
                });
        },
        getRecognition() {
            axios.post(`http://127.0.0.1:8000/api/recognition`, {
                frame: this.$store.state.frame,
            })
                .then(res => {
                    this.yData.recognition = [];
                    this.xData.recognition = [];

                    for (let x in res.data.recognition) {
                        this.yData.recognition.push(res.data.recognition[x]);
                        this.xData.recognition.push(parseInt(x) + 1);
                    }
                    this.loaded.recognition = true;
                });
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getWeight();
        this.getBulkiness();
        this.getRecognition();
    },
    mounted() {
        this.labels.weight = this.$t('charts.mass');
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