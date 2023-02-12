<template>
    <div class="drop-down">

        <p>{{ name }}</p>
        <div class="drop-down-content">
            <router-link v-for="element in hrefs" :to="`${element.path}`">{{ element.name }}</router-link>
        </div>
    </div>
</template>

<script>

export default {
    name: 'DropDownMenu',
    props: {
        name: String,
        hrefs: Array

    },
    data() {
        return {
            a: 'aaa'
        }
    }
}
</script>

<style scoped>
p {
    margin: 0
}

a {
    color: var(--text-color);
    text-decoration: none;
}


.drop-down-content {
    display: none;
    position: absolute;
    z-index: 20;
    background-color: var(--background-color);
    border: 1px solid var(--accent-color);
    border-top: 0;
    left: 50%;
    transform: translateX(calc(-50% + 0.4rem));
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
    background-color: var(--accent-color);
}

.drop-down {
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
    background-image: url('../../assets/arrow.svg');
    filter: var(--icon-filter);
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