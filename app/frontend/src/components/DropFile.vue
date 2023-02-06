<template>
    <div class="main">
        <div class="dropzone-container" @dragover="dragover" @dragleave="dragleave" @drop="drop">
            <input type="file" multiple name="file" id="fileInput" class="hidden-input" @change="onChange" ref="file"
                accept=".txt, .fasta, .pdb" />


            <label for="fileInput" class="file-label" v-if="files.length == 0">
                <div v-if="isDragging">Upuść tutaj żeby załadować.</div>
                <div v-else>Upuść Plik lub <u>Wybierz</u> żeby załadować</div>
                <img src="@/assets/upload.svg">
            </label>
            <div class="preview-container mt-4" v-if="files.length">
                <div v-for="file in files" :key="file.name" class="preview-card">
                    <div>
                        <p>
                            <img src="@/assets/file.svg">{{ file.name }} ({{ (file.size / 1000).toFixed(2) + "KB" }})

                        </p>
                    </div>
                    <div>
                        <div class="delete" type="button" @click="remove(files.indexOf(file))" title="Remove file">
                            <b>x</b>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
</template>

<script>
export default {
    name: 'DropFile',
    data() {
        return {
            isDragging: false,
            files: [],
        };
    },
    methods: {
        remove(i) {
            this.files.splice(i, 1);
        },
        onChange() {
            this.files = [...this.$refs.file.files];
        },
        dragover(e) {
            e.preventDefault();
            this.isDragging = true;
        },
        dragleave() {
            this.isDragging = false;
        },
        drop(e) {
            e.preventDefault();
            this.$refs.file.files = e.dataTransfer.files;
            this.onChange();
            this.isDragging = false;
        },
    },
};
</script>

<style scoped>
p {
    display: flex;
    align-items: center;
    justify-content: center;
    letter-spacing: 1px;
    font-size: 1.25rem;
    margin: 0;
    height: fit-content;
}

p img {
    margin-left: 0.5rem;
}

.file-label img {
    filter: var(--icon-filter);
    width: 4rem;
    height: 4rem;
}

.main {
    display: flex;
    position: relative;
    flex-grow: 1;
    align-items: center;
    height: 100%;
    width: 100%;
    justify-content: center;
    text-align: center;
}

.delete {
    top: -0.25rem;
    right: -1rem;
    cursor: pointer;
    padding: 0.25rem;
    position: absolute;

}

.dropzone-container {
    width: 100%;
    height: 100%;
}

.hidden-input {
    opacity: 0;
    overflow: hidden;
    position: absolute;
    width: 1px;
    height: 1px;
}

.file-label {
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    font-size: 20px;
    display: flex;
    cursor: pointer;
}

.file-label {
    display: flex;
    justify-content: center;
    align-items: center;
    flex-direction: column;
}

.preview-container {
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);


}

.preview-card {
    display: flex;
    position: relative;
    align-items: center;
    justify-content: space-between;

    border-bottom: rgba(0, 0, 0, 0.2) 1px solid;
    padding: 5px;

}

.preview-img {
    width: 50px;
    height: 50px;
    border-radius: 5px;
    border: 1px solid var(--accent-color-dark);
    background-color: var(--accent-color-light);
}
</style>