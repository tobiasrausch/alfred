library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
x=read.table(args[1], header=T, comment.char="$")

# Plot metrics across samples
pdf(paste0(args[1], ".pdf"))
for (s in c("ErrorRate", "DuplicateFraction", "MappedFraction", "RatioMapped2vsMapped1", "MappedSameChrFraction", "MedianReadLength", "MedianInsertSize", "MedianCoverage")) {
    print(s)
    p1=ggplot(data=x, aes_string(x="Sample", y=s))
    p1=p1 + geom_boxplot(aes(group=Library, colour=Sample), show.legend=F)
    p1=p1 + xlab("Sample") + ylab(s)
    p1=p1 + coord_flip()
    print(p1)
}
dev.off()
print(warnings())
