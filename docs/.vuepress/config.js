module.exports = {
  title: "Alfred documentation",
  description:
    "Documentation of Alfred, an app for BAM alignment statistics, feature counting and annotation",
  themeConfig: {
    nav: [
      { text: "Home", link: "/" },
      { text: "Installation", link: "/installation/" },
      { text: "Usage", link: "/cli/" },
      { text: "Web App", link: "/webapp/" }
    ],
    sidebar: ["/installation/", "/cli/", "/webapp/"]
  }
};
