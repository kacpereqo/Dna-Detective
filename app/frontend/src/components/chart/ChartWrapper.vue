<template>
    <div class="chart-wrapper" v-observe-visibility="visibilityChanged">
        <Line :word="title" />
        <div class="buttons">
            <ExportChart v-if="mounted" :yData="yData" :xData="xData" :labels="labels" :wholeNumbers="wholeNumbers"
                ref="chart" />
        </div>
        <div class="charts">
            <div id="chart">
                <Chart v-if="isVisible" :yData="yData" :xData="xData" :labels="labels" :parent="parent"
                    :wholeNumbers="wholeNumbers" ref="chart" />
            </div>
        </div>

    </div>
</template>

<script>
import ExportChart from "@/components/chart/ExportChart.vue";
import Chart from '@/components/chart/Chart.vue'
import Line from '@/components/other/Line.vue'

export default {
    name: "ChartWrapper",
    components: {
        Chart,
        ExportChart,
        Line
    },

    data() {
        return {
            isVisible: false,
            parent: null,
            mounted: false,
            title: '',
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
            type: String,
            required: false,
            default: '',
        },
        wholeNumbers: {
            type: Boolean,
            required: false,
            default: false,
        }
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
        this.title = this.$t(this.labels + '.title')
    },
    watch: {
        '$i18n.locale': function (newVal, oldVal) {
            this.title = this.$t(this.labels + '.title')
        }
    }

}

</script>

<style scoped lang="scss">
.chart-wrapper {
    width: 100%;
    height: 250px;
    margin: 2rem 0 6rem 0;
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
