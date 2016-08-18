library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
ins=read.table(args[1], header=T)

# Remove 0s
ins=ins[ins$Count!=0,]

# Add pseudo-count for log-scale
ins$InsertSize = ins$InsertSize + 1

# Generate title
count=c()
titleVal="Insert Size ("
for (o in unique(ins$Layout)) {
    s=sum(ins[ins$Layout==o,]$Count)
    count=c(count, s)
    lab=paste("#", o, ":", s, sep="")
    titleVal=paste(titleVal, lab, sep="  ")
}
countTotal=sum(count)
titleVal=paste0(titleVal, "  )")

# Plot
pdf(paste0(args[1], ".pdf"))
p1=ggplot(data=ins, aes(x=InsertSize, y=Count))
p1=p1 + geom_freqpoly(aes(colour=Layout), size=1.2, stat="identity")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle(titleVal)
p1=p1 + scale_x_log10(labels=comma)
p1=p1 + facet_wrap(Sample ~ Library)
p1
print(warnings())

