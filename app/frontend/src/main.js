import { createApp } from 'vue'
import App from './App.vue'
import router from './router'
import store from './store'
import { createMetaManager } from 'vue-meta'
import { ObserveVisibility } from 'vue-observe-visibility'
import i18n from './i18n'


const app = createApp(App);

app.directive('observe-visibility', {
    beforeMount: (el, binding, vnode) => {
        vnode.context = binding.instance;
        ObserveVisibility.bind(el, binding, vnode);
    },
    update: ObserveVisibility.update,
    unmounted: ObserveVisibility.unbind,
});

app.directive('click-outside', {
    mounted(el, binding, vnode) {
        el.clickOutsideEvent = function (event) {
            if (!(el === event.target || el.contains(event.target))) {
                binding.value(event, el);
            }
        };
        document.body.addEventListener('click', el.clickOutsideEvent);
    },
    unmounted(el) {
        document.body.removeEventListener('click', el.clickOutsideEvent);
    }
});

app.use(i18n);
app.use(createMetaManager())
app.use(store);
app.use(router);
app.mount('#app');