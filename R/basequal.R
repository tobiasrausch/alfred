library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
bq=read.table(args[1], header=T)

# Plot
png(paste0(args[1], ".png"), height=400, width=800)
p1=ggplot(data=bq, aes(x=Position, y=BaseQual))
p1=p1 + geom_line(aes(group=Library, colour=Library))
p1=p1 + xlab("Position in read") + ylab("Mean Base Quality")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Base Quality Distribution")
p1=p1 + scale_x_continuous(labels=comma)
p1=p1 + theme(legend.position="bottom", legend.direction='vertical')
p1
dev.off()
print(warnings())

