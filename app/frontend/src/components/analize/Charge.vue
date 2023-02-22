<template>
    <div class="property-wrapper">
        <h2>Ładunek</h2>
        <li> Punkt izoelektryczny {{ isoelectricpoint }}</li>
        <ChartWrapper v-if="loaded" :yData="yData" :xData="xData" :labels="labels" />
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
            xData: [],
            yData: [],
            labels: '',
            loaded: false,
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getCharge();
        this.getIsoelectricPoint();
    },
    mounted() {
        this.labels = this.$t('charts.chargeAtPH');
    },

    methods: {
        getCharge() {
            axios.get(`http://127.0.0.1:8000/api/netcharge/${this.id}?start=0&end=14&step=0.1`)
                .then(response => {
                    for (let x in response.data.netcharge) {
                        this.yData.push(response.data.netcharge[x]);
                        this.xData.push(parseFloat(x));
                        this.loaded = true;
                    }

                })
                .catch(error => {
                    console.log(error);
                })
        },
        getIsoelectricPoint() {
            axios.get(`http://127.0.0.1:8000/api/isoelectricpoint/${this.id}?scale=Rodwell`)
                .then(response => {
                    this.isoelectricpoint = response.data.isoelectricpoint;
                })
                .catch(error => {
                    console.log(error);
                })
        },
    }
}

</script>


<style ></style>