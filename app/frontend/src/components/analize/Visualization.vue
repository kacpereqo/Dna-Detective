<template>
    <div class="property-wrapper">


        <div v-if="isLoading">
            Ładowanie...
        </div>


        <Line :word="'Wzór szkieletowy'" :line="false" />

        <div class="image-container">
            <div class="fullscreen">
                Pełny Ekran<img src="@/assets/fullscreen.svg">
            </div>

            <div class="image-viewer"><img :src='visualizationSrc' alt="Wzór szkieletowy białka"></div>
        </div>
    </div>
</template>

<script>
import Line from '../other/Line.vue';
import axios from 'axios'

export default {
    name: 'Visualization',
    components: {
        Line
    },
    data() {
        return {
            visualizationSrc: '',
            isLoading: true,
        }

    },
    methods: {
        getVisualizationSrc() {
            axios.get(`http://127.0.0.1:8000/api/visualizaiton/${this.id}`,
                {
                    responseType: 'blob',
                }
            )
                .then(res => {
                    const disposition = res.headers.get('Content-Disposition');
                    var filename = disposition.split(/;(.+)/)[1].split(/=(.+)/)[1];
                    if (filename.toLowerCase().startsWith("utf-8''"))
                        filename = decodeURIComponent(filename.replace("utf-8''", ''));
                    else
                        filename = filename.replace(/['"]/g, '');
                    return res.data;
                }).then(blob => {
                    this.visualizationSrc = window.URL.createObjectURL(blob);
                    this.isLoading = false;
                });
        }

    },
    created() {
        this.id = this.$route.params.id;
        this.getVisualizationSrc();
    },

}

</script>

<style scoped>
.image-viewer img {
    filter: var(--visualization-filter);
    height: 200px;
    position: absolute;
    top: 50%;
    transform: translate(0, -50%);
}

.image-viewer {
    display: flex;
    height: 200px;
    overflow: scroll;
    overflow-y: hidden;
    scrollbar-width: thin;
    position: relative;
}

.image-container {
    margin-top: 2rem;
    position: relative;
}

.fullscreen {
    color: var(--text-color);
    border: 1px solid var(--accent-color-dark);
    width: fit-content;
    text-align: center;
    position: absolute;
    right: 0;
    padding: 0.25rem 1.75rem 0.25rem 0.5rem;
    border-radius: 0.3rem;
    cursor: pointer;
    top: -1rem;
    z-index: 10;
}

.fullscreen img {
    filter: var(--icon-filter);
    width: 1.5rem;
    height: 1.5rem;
    top: 50%;
    position: absolute;
    transform: translate(0, -50%);

}

.fullscreen:hover {
    background-color: var(--accent-color-light);
    transition: 0.2s;
}
</style>