library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
mq=read.table(args[1], header=T)
mq$RG = substr(mq$Library, 0, 12)

# Plot
pdf(paste0(args[1], ".pdf"), height=4, width=8)
p1=ggplot(data=mq, aes(x=MappingQuality, y=Fraction))
p1=p1 + geom_line(aes(group=RG, colour=RG))
p1=p1 + xlab("Mapping Quality") + ylab("Fraction of reads")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Mapping Quality Distribution")
p1=p1 + scale_x_continuous(labels=comma)
p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
p1
dev.off()
print(warnings())

