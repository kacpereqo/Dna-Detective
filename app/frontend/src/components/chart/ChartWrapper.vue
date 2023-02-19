<template>
    <div class="chart-wrapper" v-observe-visibility="visibilityChanged">
        <div class="buttons">
            <ExportChart v-if="mounted" :data="data" :labels="labels" :wholeNumbers="wholeNumbers" :xUnit="xUnit"
                ref="chart" :yUnit="yUnit" :info="' 123'" />
        </div>
        <div class="charts">

            <div id="chart">
                <Chart v-if="isVisible" :data="data" :labels="labels" :element="element" :wholeNumbers="wholeNumbers"
                    :xUnit="xUnit" ref="chart" :yUnit="yUnit" />
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
            element: null,
            mounted: false,
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
        wholeNumbers: {
            type: Boolean,
            required: false,
            default: false,
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
        this.element = this.$el.querySelector("#chart")
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
    width: calc(100vw - 216px - 3rem);
    height: 250px;
    margin-bottom: 2rem;
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
