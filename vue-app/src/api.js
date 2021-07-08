import axios from 'axios'

const jsonAxios = axios.create({
  // baseURL: 'http://localhost:5000', // required for development of vue app
  withCredentials: false,
  headers: {
    'Accept': 'application/json', // https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Accept
    'Content-Type': 'application/json' // https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Type
  }
})

export function api(config) {
  return jsonAxios.post('/api', config)
}