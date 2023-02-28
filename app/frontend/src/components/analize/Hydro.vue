<template>
    <div class="property-wrapper">
        <ChartWrapper v-if="loaded" :yData="yData" :xData="xData" :wholeNumbers="true" :labels="labels" />
        <li> Średnia hydrofobowość {{ gravy }}</li>
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
            yData: [],
            xData: [],
            labels: '',
            loaded: false,
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getHydrophobicity();
        this.getGRAVY();
    },
    mounted() {
        this.labels = this.$t('charts.hydro');
    },

    methods: {
        getHydrophobicity() {
            axios.post(`http://127.0.0.1:8000/api/hydrophobicity?scale=Kyte-Doolittle`, {
                frame: this.$store.state.frame,
            })
                .then(response => {

                    this.yData = response.data.hydrophobicity;

                    for (let i = 0; i < this.yData.length; i++) {
                        this.xData.push(i + 1);
                    }
                    this.loaded = true;
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
    }
}

</script>