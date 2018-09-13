module.exports = {
  globDirectory: 'dist/',
  globPatterns: ['**/*.{js,css,html,jpg,jpeg,png,svg,ttf,woff,woff2,eot}'],
  swDest: 'dist/sw.js',
  runtimeCaching: [
    {
      urlPattern: /^https:\/\/cdn\.plot\.ly\//,
      handler: 'staleWhileRevalidate',
      options: {
        cacheableResponse: {
          statuses: [0, 200]
        }
      }
    },
    {
      urlPattern: /\/examples\/.+\.json\.gz$/,
      handler: 'staleWhileRevalidate',
      options: {
        cacheableResponse: {
          statuses: [0, 200]
        }
      }
    },
    {
      urlPattern: /^https:\/\/fonts\.googleapis\.com/,
      handler: 'staleWhileRevalidate',
      options: {
        cacheName: 'google-fonts-stylesheets',
        cacheableResponse: {
          statuses: [0, 200]
        }
      }
    },
    {
      urlPattern: /^https:\/\/fonts\.gstatic\.com/,
      handler: 'cacheFirst',
      options: {
        cacheName: 'google-fonts-webfonts',
        cacheableResponse: {
          statuses: [0, 200]
        },
        expiration: {
          maxAgeSeconds: 60 * 60 * 24 * 365
        }
      }
    }
  ]
}
