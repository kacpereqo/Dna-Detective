<template>
    <div class="property-wrapper">
        <h2>Ładunek</h2>
        <li> Punkt izoelektryczny {{ isoelectricpoint }}</li>
        <ChartWrapper v-if="loaded" :yData="yData.charge" :xData="xData.charge" :labels="labels.charge" />
        <ChartWrapper v-if="loaded" :yData="yData.polarity" :xData="xData.polarity" :labels="labels.polarity" />
    </div>
</template>

<script>
import ChartWrapper from '@/components/chart/ChartWrapper.vue'

import axios from 'axios'

export default {
    name: 'Charge',
    components: {
        ChartWrapper,
    },
    data() {
        return {
            isoelectricpoint: '',
            xData: {},
            yData: {},
            labels: {},
            loaded: false,
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getCharge();
        this.getIsoelectricPoint();
        this.getPolarity();
    },
    mounted() {
        this.labels.charge = this.$t('charts.chargeAtPH');
    },

    methods: {
        getCharge() {
            axios.post(`http://127.0.0.1:8000/api/netcharge?start=0&end=14&step=0.1`, {
                frame: this.$store.state.frame,
            })
                .then(response => {
                    this.yData.charge = [];
                    this.xData.charge = [];
                    for (let x in response.data.netcharge) {
                        this.yData.charge.push(response.data.netcharge[x]);
                        this.xData.charge.push(parseFloat(x));
                        this.loaded = true;
                    }

                })
                .catch(error => {
                    console.log(error);
                })
        },
        getIsoelectricPoint() {
            axios.post(`http://127.0.0.1:8000/api/isoelectricpoint?scale=Rodwell`, {
                frame: this.$store.state.frame,
            })
                .then(response => {
                    this.isoelectricpoint = response.data.isoelectricpoint;
                })
                .catch(error => {
                    console.log(error);
                })
        },
        getPolarity() {
            axios.post(`http://127.0.0.1:8000/api/polarity`, {
                frame: this.$store.state.frame,
            })
                .then(response => {
                    console.log(response.data);
                    this.yData.polarity = [];
                    this.xData.polarity = [];
                    for (let x in response.data.polarity) {
                        this.yData.polarity.push(response.data.polarity[x]);
                        this.xData.polarity.push(parseInt(x));
                        this.loaded = true;
                    }

                })
                .catch(error => {
                    console.log(error);
                })
        },
    }
}

</script>


<style ></style>