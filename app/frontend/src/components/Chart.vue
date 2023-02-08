
<script>
import uPlot from 'uplot';

export default {
    name: 'Chat',
    data() {
        return {
            chart: null,
            opts: {
                scales: {
                    x: {
                        time: false,
                        auto: true,
                        range: (u, dataMin, dataMax) => [dataMin, dataMax],
                    },
                    y: {
                        auto: true,
                        range: (u, dataMin, dataMax) => [dataMin, dataMax],
                    },
                },
                title: "",
                id: "chart1",
                class: "my-chart",
                width: this.element.offsetWidth,
                height: this.element.offsetHeight,
                axes: [
                    {
                        stroke: "rgba(255, 255, 255, 0.25)", grid: {
                            show: true,
                            stroke: "rgba(255, 255, 255, 0.05)",
                            width: 2,
                            dash: [],
                        },
                        ticks: {
                            show: true,
                            stroke: "rgba(255, 255, 255, 0.25)",
                            width: 2,
                            dash: [],
                            size: 10,
                        }
                    },
                    {
                        show: true,
                        label: "Population",
                        labelSize: 30,
                        gap: 5,
                        size: 50,
                        stroke: "rgba(255, 255, 255, 0.25)",
                        grid: {
                            show: true,
                            stroke: "rgba(255, 255, 255, 0.05)",
                            width: 2,
                            dash: [],
                        },
                        ticks: {
                            show: true,
                            stroke: "rgba(255, 255, 255, 0.25)",
                            width: 2,
                            dash: [],
                            size: 10,
                        }
                    }
                ],
                series: [
                    {
                        label: "Low",
                        fill: "rgba(0, 255, 0, .2)",
                        band: true,
                    },

                    {
                        label: "Ładuenk",
                        stroke: "blue",
                        width: 2,
                        fill: "rgba(0, 0, 255,0.3)",
                    }
                ],
            },

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
        element: {
            type: Object,
            required: true,
            default: null,
        },
    },
    methods: {
        resizeChart() {
            this.chart.setSize({
                width: this.element.clientWidth,
                height: this.element.clientHeight,
            });
        },
    },
    mounted() {
        const data = [
            this.$props.labels,
            this.$props.data,
        ]


        this.chart = new uPlot(this.opts, data, this.element);
    },
    created() {
        window.addEventListener('resize', this.resizeChart);
    },
    beforeDestroy() {
        window.removeEventListener('resize', this.resizeChart);
    },
}
</script>

<style scoped>
@import 'https://unpkg.com/uplot@1.6.24/dist/uPlot.min.css';
</style>