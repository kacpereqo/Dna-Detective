import { createRouter, createWebHashHistory } from 'vue-router'
import HomeView from '../views/HomeView.vue'

const routes = [
  {
    path: '/',
    name: 'home',
    component: HomeView
  },
  {
    path: '/analize/:id',
    name: 'analize',
    component: () => import(/* webpackChunkName: "about" */ '../views/AnalizeView.vue')

  },
  {
    path: '/translations/:id',
    name: 'translations',
    component: () => import(/* webpackChunkName: "about" */ '../views/TranslationsView.vue')
  },
  {
    path: '/login',
    name: 'login',
    component: () => import(/* webpackChunkName: "about" */ '../views/LoginView.vue')
  },
  {
    path: '/register',
    name: 'register',
    component: () => import(/* webpackChunkName: "about" */ '../views/RegisterView.vue')
  },
]

const router = createRouter({
  history: createWebHashHistory(),
  routes
})

export default router
