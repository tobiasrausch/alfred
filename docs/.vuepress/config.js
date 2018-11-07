module.exports = {
  title: "Alfred documentation",
  description:
    "Documentation of Alfred, an app for BAM alignment statistics, feature counting and annotation",
  themeConfig: {
    repo: "tobiasrausch/alfred",
    nav: [
      { text: "Home", link: "/" },
      { text: "Installation", link: "/installation/" },
      { text: "Usage", link: "/cli/" },
      { text: "Web App", link: "/webapp/" },
      { text: "FAQ", link: "/faq/" }
    ],
    sidebar: ["/installation/", "/cli/", "/webapp/", "/faq/"]
  },
  plugins: {
    "@vuepress/back-to-top": true
  }
};
