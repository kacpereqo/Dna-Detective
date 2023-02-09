<template>
    <div class="chart-wrapper" v-observe-visibility="visibilityChanged">
        <div @click="saveChart">wyeksportuj</div>
        <Chart v-if="isVisible" :data="data" :labels="labels" :element="element" :wholeNumbers="wholeNumbers"
            :xUnit="xUnit" :yTitle="yTitle" ref="chart" :yUnit="yUnit" />
        <div id="chart">

        </div>
    </div>
</template>

<script>
import domtoimage from "dom-to-image-more";
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
        saveChart() {
            const canvas = document.getElementsByClassName('u-wrap')[0];
            const copy = canvas.cloneNode(true);
            domtoimage.toPng(
                canvas, { bgcolor: getComputedStyle(document.body).getPropertyValue('--background-color') }
            ).then(
                function (dataUrl) {
                    var link = document.createElement("a");
                    link.download = "chart";
                    link.href = dataUrl;
                    link.click();
                }
            )

        }
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
    margin: 3rem 0;
}


#chart {
    width: calc(100vw - 216px - 3rem);
    height: 250px;
    margin-bottom: 2rem;
}
</style>
