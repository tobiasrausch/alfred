library(gplots)

args=commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=F, comment.char="$")

pdf(paste0(args[1], ".pdf"), width=8, height=8)
textplot(t(x), show.rownames=F, show.colnames=F)
dev.off()
print(warnings())
