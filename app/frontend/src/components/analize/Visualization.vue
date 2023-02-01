<template>
    <div class="visualization-wrapper">
        <h1>Wzór szkieletowy</h1>
        <div v-if="isLoading">
            Ładowanie...
        </div>
        <img :src='visualizationSrc'>
    </div>
</template>

<script>
import axios from 'axios'

export default {
    name: 'Visualization',
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
                    console.log(filename);
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
.visualization-wrapper {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
}

img {
    width: 100%;
    height: 100%;
}
</style>