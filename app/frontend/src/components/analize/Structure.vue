<template>
    <div class="property-wrapper">
        <ChartWrapper v-if="loaded.mutability" :yData="yData.mutability" :xData="xData.mutability" :wholeNumbers="true"
            :labels="'charts.mutability'" />
        <ChartWrapper v-if="loaded.refractivity" :yData="yData.refractivity" :xData="xData.refractivity"
            :wholeNumbers="true" :labels="'charts.refractivity'" />
    </div>
</template>

<script>
import ChartWrapper from '../chart/ChartWrapper.vue';
import axios from 'axios'

export default {
    name: 'Structure',
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
                alphaHelix: false,
                refractivity: false,
            },
        }
    },
    created() {
        this.getAlphaHelix();
        this.getMutability();

    },
    methods: {
        getAlphaHelix() {
            axios.post(`http://127.0.0.1:8000/api/refractivity`, {
                frame: this.$store.state.frame,
            })
                .then(res => {
                    this.yData.refractivity = [];
                    this.xData.refractivity = [];

                    for (let x in res.data.refractivity) {
                        this.yData.refractivity.push(res.data.refractivity[x]);
                        this.xData.refractivity.push(parseInt(x) + 1);
                    }
                    this.loaded.refractivity = true;
                });
        },
        getMutability() {
            axios.post(`http://127.0.0.1:8000/api/mutability`, {
                frame: this.$store.state.frame,
            })
                .then(res => {
                    this.yData.mutability = [];
                    this.xData.mutability = [];

                    for (let x in res.data.mutability) {
                        this.yData.mutability.push(res.data.mutability[x]);
                        this.xData.mutability.push(parseInt(x) + 1);
                    }
                    this.loaded.mutability = true;
                });
        }
    }
}

</script>