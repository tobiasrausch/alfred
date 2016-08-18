library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
cov=read.table(args[1], header=T)

# Remove 0s
cov=cov[cov$Count!=0,]

# Generate title
count=c()
titleVal=paste0("Coverage Distribution")

# Plot
pdf(paste0(args[1], ".pdf"))
p1=ggplot(data=cov, aes(x=Coverage, y=Count))
p1=p1 + geom_freqpoly(aes(colour=Library), size=1.2, stat="identity")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle(titleVal)
p1=p1 + scale_x_log10(labels=comma)
p1=p1 + facet_wrap(~ Sample)
p1
print(warnings())

