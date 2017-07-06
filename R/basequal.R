library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
bq=read.table(args[1], header=T)
bq$RG = substr(bq$Library, 0, 12)

# Plot
pdf(paste0(args[1], ".pdf"), height=4, width=8)
p1=ggplot(data=bq, aes(x=Position, y=BaseQual))
p1=p1 + geom_line(aes(group=RG, colour=RG))
p1=p1 + xlab("Position in read") + ylab("Mean Base Quality")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Base Quality Distribution")
p1=p1 + scale_x_continuous(labels=comma)
p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
p1
dev.off()
print(warnings())

