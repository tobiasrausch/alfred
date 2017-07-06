library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
base=read.table(args[1], header=T)
base$RG = substr(base$Library, 0, 12)

# Plot
pdf(paste0(args[1], ".pdf"), height=4, width=8)
p1=ggplot(data=base, aes(x=Position, y=Fraction))
p1=p1 + geom_line(aes(group=Base, colour=Base))
p1=p1 + xlab("Position in read") + ylab("Base Content Fraction")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Base Content Distribution")
p1=p1 + scale_x_continuous(labels=comma)
p1=p1 + facet_wrap(~ RG) + theme(legend.position="bottom", legend.direction='horizontal')
p1
dev.off()
print(warnings())

