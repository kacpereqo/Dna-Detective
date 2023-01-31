<template>
    <div class="charge-wrapper">
        <h1>Ładunek</h1>
        <li> Punkt izoelektryczny {{ isoelectricpoint }}</li>
        <v-frappe-chart type="line" :labels="labels" :data="[{ values: charge }]" :colors="['red']" :axisOptions="{
            xIsSeries: true,
            xAxisMode: 'tick', yAxisMode: 'span'
        }" :lineOptions="{
    hideDots: 1,   //default:0
    regionFill: 1
}" />
    </div>
</template>

<script>
import { VFrappeChart } from "vue-frappe-chart"
import axios from 'axios'

export default {
    name: 'Charge',
    components: {
        VFrappeChart,
    },
    data() {
        return {
            isoelectricpoint: '',
            charge: [],
            labels: [],
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
.charge-wrapper {
    width: 80%;
}

.charge-wrapper .chart-container::after,
.chart-container::before {
    content: "Ładunek";
    font-size: 1.25rem;
    display: block;

    top: 50%;
    left: 0;
    position: absolute;
    transform: rotate(-90deg) translate(-50%, -50%);
    transform-origin: left top 0;
}

.chart-container::before {
    content: "pH";
    transform: rotate(0deg) translate(0, 90%);

    bottom: 0;
    left: 50%;

}
</style>