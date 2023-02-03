<template>
    <div class="property-wrapper">
        <h2>Ładunek</h2>
        <li> Punkt izoelektryczny {{ isoelectricpoint }}</li>
        <ChartWrapper :data="charge" :labels="labels" :lineOptions="lineOptions" />
    </div>
</template>

<script>
import ChartWrapper from '@/components/ChartWrapper.vue'

import axios from 'axios'

export default {
    name: 'Charge',
    components: {
        ChartWrapper,
    },
    data() {
        return {
            isoelectricpoint: '',
            charge: [],
            labels: [],
            lineOptions: {
                regionFill: 1,
                hideDots: 1
            },
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getCharge();
        this.getIsoelectricPoint();
    },

    methods: {
        getCharge() {
            axios.get(`http://127.0.0.1:8000/api/netcharge/${this.id}?start=0&end=14&step=0.1`)
                .then(response => {
                    for (let x in response.data.netcharge) {
                        this.charge.push(response.data.netcharge[x]);
                        this.labels.push(x);
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


<style >

</style>