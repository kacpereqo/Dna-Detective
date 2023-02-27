<template></template>

<script>
export default {
    name: "LogoutView",
    mounted() {
        this.$store.commit('logout')
        this.$router.push('/')
    },
}
</script>