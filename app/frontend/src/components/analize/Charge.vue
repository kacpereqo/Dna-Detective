<template>
    <div class="property-wrapper">
        <ChartWrapper v-if="loaded.charge" :yData="yData.charge" :xData="xData.charge" :labels="'charts.chargeAtPH'" />
        <p> {{ $t('isoelectricPoint') }} {{ isoelectricpoint }}</p>
        <ChartWrapper v-if="loaded.polarity" :yData="yData.polarity" :xData="xData.polarity" :labels="'charts.polarity'" />
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
            loaded: {
                charge: false,
                polarity: false,
            },
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getCharge();
        this.getIsoelectricPoint();
        this.getPolarity();
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
                    }
                    this.loaded.charge = true;

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
                    this.yData.polarity = [];
                    this.xData.polarity = [];
                    for (let x in response.data.polarity) {
                        this.yData.polarity.push(response.data.polarity[x]);
                        this.xData.polarity.push(parseInt(x));
                    }
                    this.loaded.polarity = true;
                })
                .catch(error => {
                    console.log(error);
                })
        },
    }
}

</script>


<style ></style>