<template>
    <div class="property-wrapper">


        <div v-if="isLoading">
            Ładowanie...
        </div>
        <Line :word="'Wzór szkieletowy'" :line="false" />
        <div class="image-container">
            <div class="image-viewer"><img :src='visualizationSrc' alt="Wzór szkieletowy białka"></div>
        </div>
        <Line :word="'Wzór'" :line="true" />
        <div class="image-container2">
            <div v-for="x in this.$store.state.frame">
                <img :src="getImgUrl(x)">
                <p> {{ x }} </p>
            </div>
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
            sequence: '',
        }

    },
    methods: {
        getVisualizationSrc() {
            axios.post(`http://127.0.0.1:8000/api/visualizaiton?`, {
                frame: this.$store.state.frame,
            },
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
        }, getImgUrl(pet) {
            var images = require.context('@/assets/compounds', false, /\.png$/)
            return images('./' + pet + ".png")
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
    justify-content: flex-start;
    align-items: center;
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

.image-container2 {
    width: calc(100vw - 16.5rem);
    display: flex;
    justify-content: flex-start;
    overflow: scroll;
    align-items: center;
}

.image-container2 img {
    filter: var(--visualization-filter);
}

.image-viewer2 {
    display: flex;
    justify-content: center;
    align-items: center;
    height: 200px;
    /* width: 100%;
    flex */
    width: 400px;
    overflow: hidden;
    overflow-x: scroll;
    /* overflow: scroll; */
    /* overflow-y: hidden; */
    scrollbar-width: thin;
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
</style>