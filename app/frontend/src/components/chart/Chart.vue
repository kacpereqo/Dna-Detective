<template></template>
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
            required: false,
            default: [],
        },
        labels: {
            type: Array,
            required: false,
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
        element: {
            type: Object,
            required: false,
            default: null,
        },
        wholeNumbers: {
            type: Boolean,
            required: false,
            default: false,
        },

    },
    methods: {
        tooltipPlugin(opts) {
            let over, bound, bLeft, bTop;

            function syncBounds(element) {
                let bbox = element.getBoundingClientRect();
                bLeft = bbox.left;
                bTop = bbox.top;
            }

            const overlay = document.createElement("div");
            overlay.id = "overlay";
            overlay.style.display = "none";
            overlay.style.position = "absolute";
            document.body.appendChild(overlay);

            overlay.style.background = "var(--background-color)";
            overlay.style.color = "(--text-color)";
            overlay.style.padding = "0.5rem";
            overlay.style.borderRadius = "0.2rem";
            overlay.style.border = "1px solid var(--accent-color-dark)";
            overlay.style.zIndex = "1000";
            overlay.style.whiteSpace = "break";


            return {
                hooks: {
                    init: u => {
                        over = u.over;

                        bound = over;

                        over.onmouseenter = () => {
                            overlay.style.display = "block";
                        };

                        over.onmouseleave = () => {
                            overlay.style.display = "none";
                            overlay.style.left = "-100vw";
                        };
                    },
                    setSize: u => {
                        syncBounds(this.element);
                    },
                    setCursor: u => {
                        const { left, top, idx } = u.cursor;
                        const x = u.data[0][idx];
                        const y = u.data[1][idx];
                        overlay.innerHTML = `<div> ${x}</div><div>${y}</div>`;

                        overlay.style.left = left + bLeft + 70 - overlay.clientWidth + "px";
                        overlay.style.top = top + bTop - overlay.clientHeight / 2 - 8 + "px";

                    }
                }
            };
        },
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
                    plugins: [
                        this.tooltipPlugin(this.element),
                    ],
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