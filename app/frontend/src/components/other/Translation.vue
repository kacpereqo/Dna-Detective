<template>
    <div class="translations-wrapper">
        <div class="header">
            <h1>Translacje</h1>
            <h2> Wybierz białko do analizy klikając na ciemno niebieskie elementy </h2>
        </div>
        <div class="translations">
            <ul>
                {{ translations.openReadingFrames }}
                <li v-for="(translation, y) in translations" class="frame-wrapper">
                    <h2>{{ translation.direction[0] }} <img src="@/assets/arrow-ltr.svg"> {{ translation.direction[1] }}
                        Ramka {{ (y % 3 +
                            1) }} </h2>
                    <div class="frames">
                        <span class="frame" v-for="(frame, i) in translation.translatedFrames" :key="i">
                            <span v-for="(f, j) in frame" :key="j">
                                <span v-if="f.open" class="open-frame" @click="toAnalize(f)">
                                    <span class="head">{{ f.first }}</span>
                                    <span class="tail" v-if="f.rest != ''">{{ f.rest }}</span>
                                </span>
                                <span v-else>{{ f.first }}{{ f.rest }}</span>
                                <span v-if="j != frame.length - 1">-</span>
                            </span>
                            <span v-if="i != translation.translatedFrames.length - 1">-</span>
                        </span>
                    </div>
                </li>
            </ul>
        </div>
    </div>
</template>

<script>
import axios from 'axios'
import { useRoute } from 'vue-router';

const getTranslation = async () => {
    const route = useRoute();
    const response = await axios.get(`http://127.0.0.1:8000/api/${route.params.id}/translate?is_reversed=true&is_forward=true`);
    return response.data;
}

export default {
    name: 'Translations',
    data() {
        return {
            id: '',
            translations: [],
        }
    },
    async setup() {
        function split(str, index) {
            const result = [str.slice(0, index), str.slice(index)];
            return result;
        }

        var translations = await getTranslation();
        for (let x in translations) {
            for (let y in translations[x].translatedFrames) {
                var temp = translations[x].translatedFrames[y]

                if (temp[0] != 'M') {
                    temp = split(temp, 1)
                    temp[0] = { first: temp[0][0], rest: temp[0].slice(1), 'open': true }

                    const index = temp[1].indexOf('M')
                    if (index != -1) {
                        temp.push({})
                        const splited = split(temp[1], index)
                        temp[1] = { first: splited[0][0], rest: splited[0].slice(1), open: false }
                        temp[2] = { first: splited[1][0], rest: splited[1].slice(1), open: true }
                    }
                    else {
                        if (temp[1].slice(1) !== "") {

                            temp[1] = { first: temp[1][0], rest: temp[1].slice(1), 'open': false }
                        }
                        else {
                            temp = temp.slice(0, 1)
                        }
                    }
                }
                else {
                    temp = [{ first: temp[0], rest: temp.slice(1), 'open': true }]
                }
                translations[x].translatedFrames[y] = temp

            }
        }

        return {
            translations,
        }
    },
    methods: {
        async toAnalize(frame) {
            axios.post(`http://127.0.0.1:8000/api/frameid`, {
                frame: frame.first + frame.rest

            }).then(
                res => {
                    this.$store.commit('setFrame', frame.first + frame.rest)
                    this.$router.push({
                        name: 'analize',
                        params: {
                            id: res.data.id
                        }
                    })
                }
            );

        }
    },
}
</script>

<style>
a {
    color: var(--text-color);
    text-decoration: none;
}

.frames,
.frame {
    word-break: break-all;
    font-family: 'Courier New', Courier, monospace;

}

.frames {
    padding: 0.5rem;
}

.head {
    cursor: pointer;
    padding: 0.1rem 0.2rem;
    border-radius: 0.2rem 0 0 0.2rem;
    color: white;
    background-color: var(--main-color-opacity);
}

.tail {
    cursor: pointer;
    word-break: break-all;
    margin-left: 0.1rem;
    padding: 0.1rem;
    border-radius: 0 0.2rem 0.2rem 0;
    background-color: var(--main-color);
}

.open-frame {
    border-radius: 0.2rem;
    word-break: break-all;
}

.translations-wrapper {
    line-height: 1.5;
    letter-spacing: 1px;
    word-break: break-all;
    margin: 2rem 0;
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 0 0.5rem;
    border: 2px solid var(--accent-color);
}

.header {
    word-break: break-all;
    display: flex;
    flex-direction: column;
    justify-content: center;
    background-color: var(--main-color);
    padding: 0.5rem;
    width: 100%;
}

.header h1 {
    color: white;
    font-weight: normal;
    font-size: 1.5rem;
    margin: 0;
}

.header h2 {
    font-weight: normal;
    color: white;
    font-size: 1rem;
    margin: 0;
}

ul {
    list-style: none;
    padding: 0;
    margin: 0;
}

li {
    list-style: none;
    padding: 0;
    margin: 0;
}

.frame-wrapper {
    word-break: break-all;
    width: 60vw;
    height: max-content;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
    margin: 1rem 0;
    border: 2px solid var(--accent-color);
    border-radius: 0.2rem;
}

.frame-wrapper h2 {
    display: flex;
    flex-direction: row;
    justify-content: flex-start;
    align-items: center;
    word-break: break-all;
    font-weight: 500;
    font-size: 1.1rem;
    padding: 0.5rem;
    margin: 0;
    margin-left: 1rem;
    width: 100%;
    position: relative;
}

.frame-wrapper h2 img {
    width: 1.5rem;
    height: 1rem;
    object-fit: cover;
}

.frame-wrapper h2::after {
    content: '';
    position: absolute;
    bottom: 0;
    left: calc(50% - 0.5rem);
    transform: translateX(-50%);
    width: calc(100% - 2rem);
    height: 2px;
    background-color: var(--accent-color);
}
</style>