library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
rl=read.table(args[1], header=T)

# Plot
pdf(paste0(args[1], ".pdf"), height=4, width=8)
p1=ggplot(data=rl, aes(x=Readlength, y=Fraction))
p1=p1 + geom_line(aes(group=Library, colour=Library))
p1=p1 + xlab("Read length") + ylab("Fraction of reads")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Read Length Distribution")
p1=p1 + scale_x_continuous(labels=comma)
p1=p1 + theme(legend.position="bottom", legend.direction='vertical')
p1
dev.off()
print(warnings())

