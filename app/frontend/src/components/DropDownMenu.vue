<template>
    <div class="drop-down">

        <p>{{ name }}</p>
        <div class="drop-down-content">
            <router-link v-for="href in hrefs" to="{{ href.href }}">{{ href.name }}</router-link>
        </div>
    </div>
</template>

<script>

export default {
    name: 'DropDownMenu',
    props: {
        name: String,
        hrefs: Array

    }
}
</script>

<style scoped>
p {
    margin: 0
}

a {
    text-decoration: none;
    color: black;
}

.drop-down-content {
    display: none;
    position: absolute;
    background-color: #f9f9f9;
    box-shadow: 0px 8px 16px 0px rgba(0, 0, 0, 0.2);
    border: 1px solid rgba(0, 0, 0, 0.25);
    border-top: 0;
    left: -0.5rem;
    padding: 0 0.4rem;
    z-index: 10;
    white-space: nowrap;
}

.drop-down-content a {
    display: block;
    width: 100%;
    padding: 0 0.4rem;
    text-align: center;
}

.drop-down-content a:hover {
    background-color: rgba(0, 0, 0, 0.25);
}

.drop-down {
    margin: 0 1rem;
    display: block;
    padding: 0 0.75rem;
    cursor: pointer;
    position: relative;
}

.drop-down:hover .drop-down-content {
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
}

.drop-down::after {
    content: '';
    position: absolute;
    background-image: url('@/assets/arrow.svg');
    background-size: contain;
    background-repeat: no-repeat;
    background-position: center;
    aspect-ratio: 1/1;
    top: 40%;
    right: -0.5rem;
    transform: translateY(-50%);
    height: 0.75rem;
    width: 0.75rem;
    animation: reverse-spin 0.25s forwards ease-in-out;

}

.drop-down:hover::after {
    animation: rotate 0.25s forwards ease-in-out;
}

@keyframes rotate {
    0% {
        transform: rotate(0deg);
    }

    100% {
        transform: rotate(180deg);
    }
}


@keyframes reverse-spin {
    100% {
        transform: rotate(0deg);
    }

    0% {
        transform: rotate(180deg);
    }
}
</style>