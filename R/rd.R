library(ggplot2)
library(scales)

chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX")

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
x = x[x$chr %in% chrs,]
x$chr = factor(x$chr, levels=chrs)

# Iterate samples
for (i in 5:ncol(x)) {
    sample = colnames(x)[i]
    print(sample)
    df = data.frame(chr=x$chr, start=x$start + (x$end - x$start) / 2, rd=log(x[,i] / median(x[,i]))/log(2))
    p1 = ggplot(data=df, aes(x=start, y=rd))
    p1 = p1 + geom_point(pch=21, size=0.5)
    p1 = p1 + xlab("Chromosome")
    p1 = p1 + ylab("Log2 median normalized read depth")
    p1 = p1 + scale_x_continuous(labels=comma)
    p1 = p1 + facet_grid(. ~ chr, scales="free_x", space="free_x")
    p1 = p1 + theme(axis.text.x = element_text(angle=45, hjust=1))
    ggsave(paste0(sample, ".wholegenome.pdf"), width=24, height=6)
    print(warnings())
}
