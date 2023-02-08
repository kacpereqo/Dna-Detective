<template>
    <div class="chart-wrapper" v-observe-visibility="visibilityChanged">
        <Chart v-if="isVisible" :data="data" :labels="labels" :element="element" ref="chart" />
        <div id="chart">

        </div>
    </div>
</template>

<script>

import Chart from '@/components/Chart.vue'

export default {
    name: "ChartWrapper",
    components: {
        Chart
    },
    data() {
        return {
            isVisible: false,
            element: null,
        }

    },
    props: {
        data: {
            type: Array,
            required: true,
            default: [],
        },
        labels: {
            type: Array,
            required: true,
            default: [],
        },
        axisOptions: {
            type: Object,
            required: false,
            default: {}
        },
        lineOptions: {
            type: Object,
            required: false,
            default: {
                hideDots: 1,
                regionFill: 1,
            },
        },
        xTitle: {
            type: String,
            required: false,
            default: "",
        },
        yTitle: {
            type: String,
            required: false,
            default: "",
        },
        xUnit: {
            type: String,
            required: false,
            default: "",
        },
        yUnit: {
            type: String,
            required: false,
            default: "",
        },
    },
    methods: {
        visibilityChanged(isVisible, entry) {
            if (this.isVisible == true) {
                this.$refs.chart.resizeChart()
            }

            if (this.isVisible == false && isVisible == true) {
                this.isVisible = true
            }
        },
    },
    mounted() {
        this.element = this.$el.querySelector("#chart")
    },

}

</script>

<style>
.chart-wrapper {
    width: 100%;
    height: 250px;
}


#chart {
    width: calc(100vw - 216px - 3rem);
    height: 250px
}
</style>
