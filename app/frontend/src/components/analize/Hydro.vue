<template>
    <div class="property-wrapper">
        <ChartWrapper v-if="loaded.hydro" :yData="yData.hydro" :xData="xData.hydro" :wholeNumbers="true"
            :labels="'charts.hydro'" />
        <p> Gravy: {{ gravy }}</p>
        <ChartWrapper v-if="loaded.omh" :yData="yData.omh" :xData="xData.omh" :wholeNumbers="true" :labels="'charts.omh'" />
    </div>
</template>

<script>
import ChartWrapper from '../chart/ChartWrapper.vue';
import axios from 'axios'

export default {
    name: 'Hydro',
    components: {
        ChartWrapper
    },
    data() {
        return {
            gravy: '',
            yData: {},
            xData: {},
            labels: {},
            loaded: {
                hydro: false,
                omh: false
            },
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getHydrophobicity();
        this.getGRAVY();
        this.getOMH();
    },

    methods: {
        getHydrophobicity() {
            axios.post(`http://127.0.0.1:8000/api/hydrophobicity?scale=Kyte-Doolittle`, {
                frame: this.$store.state.frame,
            })
                .then(response => {

                    this.xData.hydro = []
                    this.yData.hydro = response.data.hydrophobicity;

                    for (let i = 1; i < this.yData.hydro.length + 1; i++) {
                        this.xData.hydro.push(i + 1);
                    }
                    this.loaded.hydro = true;
                })
        },
        getGRAVY() {
            axios.post(`http://127.0.0.1:8000/api/avghydrophobicity?scale=Kyte-Doolittle`, {
                frame: this.$store.state.frame,
            })
                .then(response => {
                    this.gravy = response.data.hydrophobicity;
                })
                .catch(error => {
                    console.log(error);
                })
        },
        getOMH() {
            axios.post(`http://127.0.0.1:8000/api/omh`, {
                frame: this.$store.state.frame,
            })
                .then(res => {
                    this.yData.omh = [];
                    this.xData.omh = [];

                    for (let x in res.data.omh) {
                        this.yData.omh.push(res.data.omh[x]);
                        this.xData.omh.push(parseInt(x) + 1);
                    }
                    this.loaded.omh = true;
                });
        }
    }
}

</script>