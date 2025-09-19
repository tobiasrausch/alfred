import { defineUserConfig } from '@vuepress/cli'
import { viteBundler } from '@vuepress/bundler-vite'
import { defaultTheme } from '@vuepress/theme-default'

export default defineUserConfig({
  title: "Alfred documentation",
  description: "Documentation of Alfred, an app for BAM alignment statistics, feature counting and annotation",
  base: "/docs/alfred/",
  theme: defaultTheme({
    repo: "tobiasrausch/alfred",
    navbar: [
      { text: "Home", link: "/" },
      { text: "Installation", link: "/installation/" },
      { text: "Usage", link: "/cli/" },
      { text: "Web App", link: "/webapp/" },
      { text: "FAQ", link: "/faq/" }
    ],
    sidebar: [
      { text: "Installation", link: "/installation/" },
      { text: "Usage", link: "/cli/" },	
      { text: "Web App", link: "/webapp/" },
      { text: "FAQ", link: "/faq/" }
    ]
  }),
  bundler: viteBundler({})
})
