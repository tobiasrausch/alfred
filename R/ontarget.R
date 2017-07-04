library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
cov=read.table(args[1], header=T)
cov$Extension=as.numeric(cov$Extension)
cov$OnTarget=as.numeric(cov$OnTarget)

# Plot
png(paste0(args[1], ".png"), height=400, width=800)
p1=ggplot(data=cov, aes(x=Extension, y=OnTarget))
p1=p1 + geom_line(aes(group=Library, colour=Library))
p1=p1 + scale_y_continuous(labels=comma, limits=c(0,1)) + ggtitle("On-target rate")
p1=p1 + scale_x_continuous(labels=comma) + xlab("Left/Right Extension of target region") + ylab("Fraction of reads on-target")
p1=p1 + facet_wrap(~ Sample) + theme(legend.position="bottom", legend.direction='vertical')
p1
dev.off()
print(warnings())

