<template>

    <div class="chart-wrapper" v-observe-visibility="visibilityChanged">
        <v-frappe-chart type="line" :labels="labels" :data="[{ values: data }]" :colors="['red']" :axisOptions="{
            xIsSeries: true,
            xAxisMode: 'tick', yAxisMode: 'span'
        }" :lineOptions="lineOptions" :width="1" v-if="isVisible" />
    </div>
</template>

<script>
import { VFrappeChart } from "vue-frappe-chart"

export default {
    name: 'ChartWrapper',
    components: {
        VFrappeChart,
    },
    data() {
        return {
            isVisible: false,
        }
    },
    props: {
        data: {
            type: Array,
            required: true
        },
        labels: {
            type: Array,
            required: true
        },
        axisOptions: {
            type: Object,
            required: false
        },
        lineOptions: {
            type: Object,
            required: false
        },
        xTittle: {
            type: String,
            required: false
        },
        yTittle: {
            type: String,
            required: false
        },
    },
    methods: {
        visibilityChanged(isVisible, entry) {
            if (this.isVisible == false && isVisible == true) {
                this.isVisible = isVisible
            }
            console.log("123", isVisible, entry)
        },
    },
}

</script>

<style>
.chart-wrapper {
    width: calc(100vw - 216px - 3rem);
}
</style>