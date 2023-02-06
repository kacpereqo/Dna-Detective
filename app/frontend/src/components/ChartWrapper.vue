<template>

    <div class="chart-wrapper" v-observe-visibility="visibilityChanged" :style="{'--xTitle': `'${xTitle}'`, '--yTitle'
    :`'${yTitle}'`}">
        <v-frappe-chart type="line" :labels="labels" :data="[{ values: data }]" :colors="['blue']" :axisOptions="{
            xIsSeries: true,
            xAxisMode: 'tick', yAxisMode: 'span'
        }" :lineOptions="lineOptions" :width="1" v-if="isVisible" :TooltipOptions="{
    formatTooltipX: d => (d + ' ' + ''),
    formatTooltipY: d => (d + ' ' + '')
}" />
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
            if (this.isVisible == false && isVisible == true) {
                this.isVisible = isVisible
            }

        },
    },
}

</script>

<style>

</style>