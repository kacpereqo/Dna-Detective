<template>

    <div class="chart-wrapper">

        <v-frappe-chart type="line" :labels="labels" :data="[{ values: data }]" :colors="['red']"
            :axisOptions="axisOptions" :lineOptions="lineOptions" :width="1" />
    </div>
</template>

<script>
import { VFrappeChart } from "vue-frappe-chart"

export default {
    name: 'ChartWrapper',
    components: {
        VFrappeChart,
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
    },

}

</script>

<style>
.chart-wrapper,
.chart {
    width: 100%;
}
</style>