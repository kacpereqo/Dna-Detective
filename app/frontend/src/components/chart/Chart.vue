<template></template>
<script>
import uPlot from 'uplot';

export default {
    name: 'Chat',
    data() {
        return this.initState();
    },
    props: {
        xData: {
            type: Array,
            required: false,
            default: [],
        },
        yData: {
            type: Array,
            required: false,
            default: [],

        },
        labels: {
            type: String,
            required: false,
            default: "",
        },
        parent: {
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
                        syncBounds(this.parent);
                        overlay.display = "none";
                    },
                    setCursor: u => {
                        const { left, top, idx } = u.cursor;
                        const x = u.data[0][idx];
                        overlay.style.left = left + bLeft - overlay.clientWidth + 60 + "px";
                        overlay.style.top = top + bTop - overlay.clientHeight + 20 + "px";

                        const y = u.data[1][idx];
                        const aa = this.wholeNumbers ? `<b>(${this.$store.state.frame[idx]})<b>` : ""
                        console.log(this.wholeNumbers);

                        overlay.innerHTML = `
                        <ul>
                        <li class="x"> <b> ${x} ${this.$t(this.labels + '.xUnit')} </b> ${aa} </li>
                        <li class="y"> <div style="--color:${u.series[1]._stroke}" class="square" ></div>${y} ${this.$t(this.labels + '.yUnit')}</li>
                    </ul>`;

                    }
                }
            };
        },
        resizeChart() {
            this.chart.setSize({
                width: this.parent.clientWidth,
                height: this.parent.clientHeight,
            });
        },
        reInit() {
            const data = [
                this.xData,
                this.yData,
            ]
            Object.assign(this.$data, this.initState());
            this.parent.innerHTML = "";
            this.chart = new uPlot(this.opts, data, this.parent);
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
                        this.tooltipPlugin(),
                    ],
                    width: this.parent.offsetWidth,
                    height: this.parent.offsetHeight,
                    axes: [
                        {
                            label: this.$t(this.labels + '.xAxis'),
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
                            label: this.$t(this.labels + '.yAxis'),
                            yDataize: 30,
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
                        {},
                        {
                            stroke: "blue",
                            width: 2,
                            fill: root.getPropertyValue('--chart-color'),
                        }

                    ],
                    legend: {
                        show: false,
                    }
                },

            }
        }
    },
    mounted() {
        const data = [
            this.xData,
            this.yData,
        ]
        this.chart = new uPlot(this.opts, data, this.parent);
    },
    created() {
        this.unwatch = this.$store.watch(
            (state) => state.theme,
            () => {
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

<style >
@import 'https://unpkg.com/uplot@1.6.24/dist/uPlot.min.css';

#overlay {
    position: absolute;
    background: var(--background-color);
    color: var(--text-color);
    padding: 0.5rem;
    border-radius: 0.2rem;
    border: 1px solid var(--accent-color-dark);
    z-index: 1000;
    min-width: 3rem;
}

#overlay::after {
    content: "";
    position: absolute;
    top: 50%;
    left: 100%;
    border-width: 5px;
    border-style: solid;
    border-color: transparent transparent transparent var(--text-color);
}

#overlay ul {
    list-style: none;
    padding: 0;
    margin: 0;
}

#overlay li {
    display: flex;
    align-items: center;
    margin: 0.2rem 0;
}

.square {
    width: 0.75rem;
    height: 0.75rem;
    border: 2px solid var(--accent-color);
    border-radius: 0.1rem;
    margin-right: 0.5rem;
    background-color: var(--color);
}

.x {
    font-weight: bold;
    display: flex;
    justify-content: space-between;
    border-bottom: 2px solid var(--accent-color-dark);
}

.y {
    font-weight: normal;
}
</style>