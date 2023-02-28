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

                        var placement = function () { "use strict"; var e = { size: ["height", "width"], clientSize: ["clientHeight", "clientWidth"], offsetSize: ["offsetHeight", "offsetWidth"], maxSize: ["maxHeight", "maxWidth"], before: ["top", "left"], marginBefore: ["marginTop", "marginLeft"], after: ["bottom", "right"], marginAfter: ["marginBottom", "marginRight"], scrollOffset: ["pageYOffset", "pageXOffset"] }; function t(e) { return { top: e.top, bottom: e.bottom, left: e.left, right: e.right } } return function (o, r, f, a, i) { void 0 === f && (f = "bottom"), void 0 === a && (a = "center"), void 0 === i && (i = {}), (r instanceof Element || r instanceof Range) && (r = t(r.getBoundingClientRect())); var n = Object.assign({ top: r.bottom, bottom: r.top, left: r.right, right: r.left }, r), s = { top: 0, left: 0, bottom: window.innerHeight, right: window.innerWidth }; i.bound && ((i.bound instanceof Element || i.bound instanceof Range) && (i.bound = t(i.bound.getBoundingClientRect())), Object.assign(s, i.bound)); var l = getComputedStyle(o), m = {}, b = {}; for (var g in e) m[g] = e[g]["top" === f || "bottom" === f ? 0 : 1], b[g] = e[g]["top" === f || "bottom" === f ? 1 : 0]; o.style.position = "absolute", o.style.maxWidth = "", o.style.maxHeight = ""; var d = parseInt(l[b.marginBefore]), c = parseInt(l[b.marginAfter]), u = d + c, p = s[b.after] - s[b.before] - u, h = parseInt(l[b.maxSize]); (!h || p < h) && (o.style[b.maxSize] = p + "px"); var x = parseInt(l[m.marginBefore]) + parseInt(l[m.marginAfter]), y = n[m.before] - s[m.before] - x, z = s[m.after] - n[m.after] - x; (f === m.before && o[m.offsetSize] > y || f === m.after && o[m.offsetSize] > z) && (f = y > z ? m.before : m.after); var S = f === m.before ? y : z, v = parseInt(l[m.maxSize]); (!v || S < v) && (o.style[m.maxSize] = S + "px"); var w = window[m.scrollOffset], O = function (e) { return Math.max(s[m.before], Math.min(e, s[m.after] - o[m.offsetSize] - x)) }; f === m.before ? (o.style[m.before] = w + O(n[m.before] - o[m.offsetSize] - x) + "px", o.style[m.after] = "auto") : (o.style[m.before] = w + O(n[m.after]) + "px", o.style[m.after] = "auto"); var B = window[b.scrollOffset], I = function (e) { return Math.max(s[b.before], Math.min(e, s[b.after] - o[b.offsetSize] - u)) }; switch (a) { case "start": o.style[b.before] = B + I(n[b.before] - d) + "px", o.style[b.after] = "auto"; break; case "end": o.style[b.before] = "auto", o.style[b.after] = B + I(document.documentElement[b.clientSize] - n[b.after] - c) + "px"; break; default: var H = n[b.after] - n[b.before]; o.style[b.before] = B + I(n[b.before] + H / 2 - o[b.offsetSize] / 2 - d) + "px", o.style[b.after] = "auto" }o.dataset.side = f, o.dataset.align = a } }();

                        const y = u.data[1][idx];
                        const aa = this.wholeNumbers ? `<b>(${this.$store.state.frame[idx]})<b>` : ""

                        overlay.innerHTML = `
                        <ul>
                        <li class="x"> <b> ${x} ${this.$t(this.labels + '.xUnit')} </b> ${aa} </li>
                        <li class="y"> <div style="--color:${u.series[1]._stroke}" class="square" ></div>${y} ${this.$t(this.labels + '.yUnit')}</li>
                    </ul>`;


                        const anchor = { left: left + bLeft, top: top + bTop };
                        placement(overlay, anchor, "right", "start", { bound });
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
    watch: {
        '$i18n.locale': function (newVal, oldVal) {
            console.log("locale changed");
            this.reInit();
        }
    }
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
    pointer-events: none;
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