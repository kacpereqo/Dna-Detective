
<script>
import uPlot from 'uplot';

export default {
    name: 'Chat',
    data() {
        return this.initState();
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
        wholeNumbers: {
            type: Boolean,
            required: false,
            default: false,
        },

    },
    methods: {
        resizeChart() {
            this.chart.setSize({
                width: this.element.clientWidth,
                height: this.element.clientHeight,
            });
        },
        reInit() {
            const data = [
                this.$props.labels,
                this.$props.data,
            ]
            Object.assign(this.$data, this.initState());
            this.element.innerHTML = "";
            this.chart = new uPlot(this.opts, data, this.element);
        },
        initState() {
            const root = getComputedStyle(document.body);

            return {
                chart: null,
                opts: {
                    scales: {
                        x: {
                            distr: this.wholeNumbers ? 2 : 1,
                            time: false,
                            auto: false,
                        },
                        y: {
                            auto: true,

                        },
                    },
                    title: "",
                    id: "chart1",
                    class: "my-chart",
                    width: this.element.offsetWidth,
                    height: this.element.offsetHeight,
                    axes: [
                        {
                            label: this.xUnit,
                            labelFont: "16px Arial ",
                            stroke: root.getPropertyValue('--text-color'), grid: {
                                show: true,
                                stroke: root.getPropertyValue('--accent-color-light'),
                                width: 2,
                                dash: [],
                            },
                            ticks: {
                                show: true,
                                stroke: root.getPropertyValue('--accent-color-dark'),
                                width: 2,
                                dash: [],
                                size: 10,
                            }
                        },
                        {
                            labelFont: "16px Arial ",
                            show: true,
                            label: this.yUnit,
                            labelSize: 30,
                            gap: 5,
                            size: 50,
                            stroke: root.getPropertyValue('--text-color'),
                            grid: {
                                show: true,
                                stroke: root.getPropertyValue('--accent-color-light'),
                                width: 2,
                                dash: [],
                            },
                            ticks: {
                                show: true,
                                stroke: root.getPropertyValue('--accent-color-dark'),
                                width: 2,
                                dash: [],
                                size: 10,
                            }
                        }
                    ],
                    series: [
                        {
                            label: this.xUnit,


                        },

                        {
                            label: this.yUnit,
                            stroke: "blue",
                            width: 2,
                            fill: root.getPropertyValue('--chart-color'),
                        }
                    ],
                },

            }
        }
    },
    mounted() {
        const data = [
            this.$props.labels,
            this.$props.data,
        ]

        this.chart = new uPlot(this.opts, data, this.element);
    },
    created() {
        this.unwatch = this.$store.watch(
            (state) => state.theme,
            (theme) => {
                this.reInit();
            }
        );
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