import { createRouter, createWebHashHistory } from 'vue-router'
import HomeView from '../views/HomeView.vue'

const routes = [
  {
    path: '/',
    name: 'home',
    component: HomeView
  },
  {
    path: '/analize',
    name: 'analize',
    component: () => import(/* webpackChunkName: "analize" */ '../views/AnalizeView.vue')
  },
]

const router = createRouter({
  history: createWebHashHistory(),
  routes
})

export default router
