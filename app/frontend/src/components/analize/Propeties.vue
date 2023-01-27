<template>
    <div class="propeties-wrapper">
        <p>{{ weight }} kilodaltons</p>

    </div>
</template>

<script>
import axios from 'axios'

export default {
    name: 'Propeties',
    data() {
        return {
            weight: '',
        }
    },
    methods: {
        getWeight() {
            axios.get(`http://127.0.0.1:8000/api/weight/${this.id}`)
                .then(res => {
                    this.weight = res.data.weight;
                })
        }
    },
    created() {
        this.id = this.$route.params.id;
        this.getWeight();
    },

}

</script>

<style scoped>
.propeties-wrapper {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
}
</style>