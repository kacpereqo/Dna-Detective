<template>
    <div class="chart-wrapper" v-observe-visibility="visibilityChanged">
        <div class="buttons">
            <ExportChart v-if="mounted" :yData="yData" :xData="xData" :labels="labels" ref="chart" />
        </div>
        <div class="charts">

            <div id="chart">
                <Chart v-if="isVisible" :yData="yData" :xData="xData" :labels="labels" :parent="parent" ref="chart" />
            </div>
        </div>

    </div>
</template>

<script>
import ExportChart from "@/components/chart/ExportChart.vue";
import Chart from '@/components/chart/Chart.vue'

export default {
    name: "ChartWrapper",
    components: {
        Chart,
        ExportChart
    },

    data() {
        return {
            isVisible: false,
            parent: null,
            mounted: false,
        }
    },
    props: {
        xData: {
            type: Array,
            required: true,
            default: [],
        },
        yData: {
            type: Array,
            required: true,
            default: [],
        },
        labels: {
            type: Object,
            required: false,
            default: {},
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
        reInit() {
            this.$refs.chart.reInit()
        },

    },
    mounted() {
        this.mounted = true
        this.parent = this.$el.querySelector("#chart")
    },

}

</script>

<style scoped lang="scss">
.chart-wrapper {
    width: 100%;
    height: 250px;
    margin: 3rem 0;
}


#chart {
    width: calc(100vw - 16.5rem);
    height: 250px;
    margin-bottom: 2rem;
    position: relative;
}

@media screen and (max-width: 960px) {
    #chart {
        width: calc(100vw - 3rem);
    }
}



.buttons {
    display: flex;
    justify-content: flex-end;
    margin-bottom: 0.25rem;
    margin-right: 1rem;
}
</style>
