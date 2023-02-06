<template>
    <div class="property-wrapper">
        <h2>Hydrofobowość</h2>
        <li> Średnia hydrofobowość {{ gravy }}</li>
        <ChartWrapper :labels="labels" :data="hydrophobicity" :xUnit="' Aminokwas'" :yUnit="''"
            :xTitle="'Pozycja Aminokwasu'" yTitle="Hydrofobowośćaaaa" />
    </div>
</template>

<script>
import ChartWrapper from '../ChartWrapper.vue';
import axios from 'axios'

export default {
    name: 'Hydro',
    components: {
        ChartWrapper
    },
    data() {
        return {
            gravy: '',
            hydrophobicity: [],
            labels: [],
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getHydrophobicity();
        this.getGRAVY();
    },

    methods: {
        getHydrophobicity() {
            axios.get(`http://127.0.0.1:8000/api/hydrophobicity/${this.id}?scale=Kyte-Doolittle`)
                .then(response => {
                    this.hydrophobicity = response.data.hydrophobicity;
                    for (let i = 2; i < this.hydrophobicity.length + 2; i++) {
                        this.labels.push(i.toString());
                    }

                })
                .catch(error => {
                    console.log(error);
                })
        },
        getGRAVY() {
            axios.get(`http://127.0.0.1:8000/api/avghydrophobicity/${this.id}?scale=Kyte-Doolittle`)
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