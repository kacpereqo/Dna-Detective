<template>
    <div class="hydro-wrapper">
        <h1>Hydrofobowość</h1>
        <li> Średnia hydrofobowość {{ gravy }}</li>
        <v-frappe-chart type="line" :labels="labels" :data="[{ values: hydrophobicity }]" :colors="['red']" />
    </div>
</template>

<script>
import { VFrappeChart } from "vue-frappe-chart"
import axios from 'axios'

export default {
    name: 'Hydro',
    components: {
        VFrappeChart,
    },
    data() {
        return {
            gravy: '',
            hydrophobicity: '',
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

<style scoped>
.hydro-wrapper {
    width: 70%;
}
</style>