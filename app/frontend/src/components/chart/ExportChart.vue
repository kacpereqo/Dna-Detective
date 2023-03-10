<template>
    <div class="buttons">
        <div class="info-button button" @click="showInfo = !showInfo" v-click-out-side="() => { showInfo = false }">
            <div> {{ $t('info') }} </div>
            <img src="@/assets/error.svg">
            <div class="info" v-if="showInfo">
                <p>{{ $t(labels + '.info') }}</p>
            </div>

        </div>

        <!-- Export buttons -->

        <div class="export-button ">
            <div class="export-button-primary button" @click="() => { showPopup = true; getSrc() }">{{ $t('export.title') }}
            </div>
            <div class="export-button-dropdown" v-click-out-side="() => { showOptions = false }"
                :style="{ 'border-radius': showOptions ? '0 0.3rem 0 0' : '0 0.3rem 0.3rem 0' }"
                @click="showOptions = !showOptions">
                <div class="arrow-container">
                    <div class="arrow" :class="{ 'animate': showOptions, 'reverse-animate': !showOptions }"> </div>
                </div>
                <div class="drop-content" v-show="showOptions">
                    <ul>
                        <li @click="() => { showPopup = true; getSrc() }">
                            Png
                        </li>
                        <li @click="saveCsv">
                            Csv
                        </li>
                        <li @click="saveJson">
                            Json
                        </li>
                    </ul>
                </div>
            </div>
        </div>
    </div>
    <!-- Popup Background  -->
    <div v-if="showPopup" @click="() => { showPopup = false; width = 800; height = 400; resizeChart() }"
        class="export-popup-background">
    </div>

    <!-- Popup  -->
    <div v-if="showPopup" class="export-popup">
        <div class="popup-close" @click="{ showPopup = false; width = 800; height = 400; resizeChart() }">X</div>
        <div class="container">

            <div class="preview-container">
                <p>{{ $t('export.preview') }}</p>
                <div class="preview">
                    <img :src="src" alt="chart">
                </div>
            </div>
            <div class="sidebar">
                <div>
                    <p>{{ $t('export.width') }}:</p> <input type="text" min="1" v-model="width"
                        v-on:change="() => { resizeChart(); getSrc() }">
                </div>

                <div>
                    <p>{{ $t('export.height') }}</p> <input type="text" min="1" v-model="height"
                        v-on:change="() => { resizeChart(); getSrc() }">
                </div>

                <div class="export" @click="saveChart()">
                    <p>{{ $t('export.title') }}</p>
                </div>
            </div>
        </div>
    </div>

    <!-- Export Chart -->
    <div>
        <div id="export-chart" ref="chartContainer">
            <Chart v-if="mounted" :xData="xData" :yData="yData" :parent="parent" :labels="labels" ref="chart" />
        </div>
    </div>
</template>

<script>
import clickOutSide from "@mahdikhashan/vue3-click-outside";
import domtoimage from "dom-to-image-more";
import Chart from "@/components/chart/Chart.vue";

export default {
    name: 'ExportChart',
    components: {
        Chart,
    },
    directives: {
        clickOutSide,
    },
    data() {
        return {
            width: 800,
            height: 400,
            mounted: false,
            showPopup: false,
            showInfo: false,
            showOptions: false,
            parent: this,
            src: null,
        }
    },
    props: {
        yData: {
            type: Array,
            required: true,
            default: [],
        },
        xData: {
            type: Array,
            required: true,
            default: [],
        },
        labels: {
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
        getSrc() {
            this.resizeChart();
            const chart = this.parent.getElementsByClassName("u-wrap")[0];
            domtoimage.toPng(
                chart, { bgcolor: getComputedStyle(document.body).getPropertyValue('--background-color') }
            ).then(
                (dataUrl) => {
                    this.src = dataUrl;
                }
            )
        },

        saveChart() {
            const chart = this.parent.getElementsByClassName("u-wrap")[0];
            domtoimage.toPng(
                chart, { bgcolor: getComputedStyle(document.body).getPropertyValue('--background-color') }
            ).then(
                (dataUrl) => {
                    const link = document.createElement("a");
                    link.download = "chart";
                    link.href = dataUrl;
                    link.click();
                }
            )
        },
        resizeChart() {
            this.parent.style.width = this.width + "px";
            this.parent.style.height = this.height + "px";
            this.$refs.chart.resizeChart();
        },
        saveCsv() {
        },
        saveCsv() {

            let csv = 'label,data\n';
            this.yData.forEach((row, i) => {
                csv += this.xData[i] + ',' + row + '\n';
            });
            const link = document.createElement("a");
            link.download = "data.csv";
            link.href = "data:text/csv;charset=utf-8," + encodeURIComponent(csv);
            link.click();
        },
        saveJson() {
            const link = document.createElement("a");
            link.download = "data.json";
            link.href = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify({
                data: this.data,
                labels: this.labels,
            }));
            link.click();
        },
    },
    mounted() {
        this.parent = this.$refs.chartContainer;
        this.mounted = true;
    }

}
</script>

<style scoped>
.export-button {
    display: flex;
}

.info {
    display: flex;
    align-items: center;
    justify-content: center;
    position: absolute;
    background: var(--background-color);
    bottom: calc(100% + 0.75rem);
    left: 50%;
    transform: translateX(-50%);
    z-index: 100;
}


.info p {
    width: max-content;
    max-width: 24rem;
    overflow: hidden;
    max-height: 9rem;
    margin: 0;
    padding: 0.25rem;
    border: 1px solid var(--accent-color);
    text-align: justify;
    font-size: 0.9rem;
    border-radius: 0.25rem;
}

.info::after {
    content: "";
    position: absolute;
    top: calc(100%);
    left: 50%;
    transform: translateX(-50%);
    border: 0.5rem solid transparent;
    border-top-color: var(--accent-color);
}

.buttons {
    display: flex;
    flex-direction: row;
    align-items: center;
    margin: 0.5rem 0;
}

.info-button {
    position: relative;
    display: flex;
    align-items: center;
    justify-content: space-between;
    cursor: pointer;
    margin: 0 0.1rem;
    border: 1px solid var(--accent-color);
    padding: 0.25rem 0.5rem;
    border-radius: 0.3rem;
}

.info-button img {
    margin-left: 0.25rem;
    width: 1.25rem;
    height: 1.25rem;
    filter: var(--icon-filter);
}

.button {
    margin-left: 0.5rem;
}

.info-button:hover {
    transition: 0.2s;
    background-color: var(--accent-color-light);
}

.export-button-dropdown,
.button {
    height: 1.5rem;
    cursor: pointer;

    border: 1px solid var(--accent-color);
}

.export-button-primary:hover {
    transition: 0.2s;
    background-color: var(--accent-color-light);
}

.export-button-primary {
    margin-right: 0.25rem;
    padding: 0.25rem 0.5rem;
    border-radius: 0.3rem 0 0 0.3rem;
}

.export-button-dropdown {
    position: relative;
    width: 1.5rem;
    height: 2rem;
    border-radius: 0 0.3rem 0.3rem 0;
}


.arrow-container:hover {
    transition: 0.2s;
    background-color: var(--accent-color-light);
}

.arrow-container {
    display: flex;
    justify-content: center;
    align-items: center;
    width: 100%;
    height: 100%;
}

.arrow {
    position: absolute;
    width: 0.75rem;
    height: 0.75rem;
    background-size: 0.75rem;
    background-image: url('@/assets/arrow.svg');
    background-repeat: no-repeat;
    background-position: center;
    filter: var(--icon-filter);
    animation: var(--animation);
    transform-origin: center;
}

.export-button ul {
    list-style: none;
    padding: 0;
    margin: 0;
    right: -1px;
    position: absolute;
    top: calc(100% + 3px);
    background-color: var(--background-color);
    border: 1px solid var(--accent-color);
    z-index: 1;
}


.export-button li {
    padding: 0.45rem 0.65rem;
    cursor: pointer;
    pointer-events: auto;
}

.export-button li:hover {
    transition: 0.3s;
    background-color: var(--accent-color-light);
}

.drop-content::after {
    content: "";
    position: absolute;
    width: 1.5rem;
    height: 4px;
    z-index: 5;
    border-right: 1px solid var(--accent-color);
    border-left: 1px solid var(--accent-color);
    left: -1px;
    background-color: var(--background-color);
}

@keyframes rotate {
    0% {
        transform: rotate(0deg);
    }

    100% {
        transform: rotate(180deg);
    }
}

@keyframes reverseRotate {
    0% {
        transform: rotate(180deg);
    }

    100% {
        transform: rotate(0deg);
    }
}

.animate {
    animation: rotate 0.2s linear forwards;
}

.reverse-animate {
    animation: reverseRotate 0.2s linear forwards;
}

.export-popup-background {
    position: fixed;
    top: 0;
    left: 0;
    background-color: rgba(0, 0, 0, 0.5);
    width: 100vw;
    height: 100vh;
    z-index: 10;
}

.export-popup {
    position: fixed;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    width: 90%;
    height: 80%;
    background-color: var(--background-color);
    border: 1px solid var(--accent-color);
    border-radius: 0.3rem;
    z-index: 100;
}

.popup-close {
    position: absolute;
    z-index: 1000;
    top: 0rem;
    right: 0.5rem;
    cursor: pointer;
    font-weight: bold;
    font-size: 1.5rem;
    color: var(--text-color);
    padding: 0.5rem;
}

.popup-close:hover {
    transition: 0.2s;
    color: var(--accent-color-dark);
}

.container {
    position: relative;
    width: 100%;
    height: 100%;
    display: flex;
}

.sidebar {
    width: 196px;
    height: 100%;
    background-color: var(--background-color);
    border-left: 1px solid var(--accent-color);
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
}

.sidebar div {
    margin: 0.5rem 0;
    display: flex;
}

#export-chart {
    position: absolute;
    width: 800px;
    height: 400px;
    z-index: -1;
}

.sidebar input {
    margin: 0 0.5rem;
    background-color: transparent;
    border: 1px solid var(--accent-color);
    border-radius: 0.3rem;
    width: 5rem;
    padding: 0.25rem;
}

.sidebar p {
    margin: 0;
}

.preview {
    display: flex;
    justify-content: center;
    align-items: center;
    position: relative;
    width: 100%;
    height: 100%;
    overflow: hidden;
}

.preview-container {
    display: flex;

    justify-content: flex-start;
    align-items: center;
    flex-direction: column;
    flex: 1;
}

.preview-container p {
    font-size: 1.25rem;
}

.preview img {
    max-width: 100%;
    max-height: 100%;
    object-fit: contain;
}

.export {
    font-size: 1.25rem;
    position: relative;
    border: 1px solid var(--accent-color);
    border-radius: 0.3rem;
    padding: 0.5rem 2.25rem 0.5rem 1rem;
    cursor: pointer;
    margin: 0.5rem 0;
}

.export:hover {
    transition: 0.2s;
    background-color: var(--accent-color-light);
}

.export p {
    margin: 0;
}

.export::after {
    content: "";
    position: absolute;
    right: 0.5rem;
    width: 1.5rem;
    height: 1.5rem;
    top: 50%;
    transform: translateY(-50%);
    background-size: 1.5rem;
    background-image: url('@/assets/export.svg');
    background-repeat: no-repeat;
    background-position: center;
    filter: var(--icon-filter);
}
</style>