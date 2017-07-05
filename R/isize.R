library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
ins=read.table(args[1], header=T)
ins$InsertSize=as.numeric(ins$InsertSize)
ins$Count=as.numeric(ins$Count)

# Remove insert size > maxvalue
tc=sum(as.numeric(ins$Count))
lastVal=max(ins$InsertSize)
sub=ins[ins$InsertSize < lastVal & ins$Count>0,]$InsertSize
upperbound=max(sub)
lowerbound=min(sub)
gr=sum(as.numeric(ins[ins$InsertSize>upperbound,]$Count))
infoMax = paste0("Insert size > ", upperbound, " (", round(100 * gr / tc, digits=2), "%)")
gr=sum(as.numeric(ins[ins$InsertSize<lowerbound,]$Count))
infoMin = paste0("Insert size < ", lowerbound, " (", round(100 * gr / tc, digits=2), "%)")
print(infoMax)
print(infoMin)
ins=ins[ins$InsertSize <= upperbound,]
ins=ins[ins$InsertSize >= lowerbound,]

# Plot
png(paste0(args[1], ".png"), height=400, width=800)
p1=ggplot(data=ins, aes(x=InsertSize, y=Count))
p1=p1 + geom_line(aes(group=Layout, colour=Layout))
p1=p1 + scale_y_continuous(labels=comma) + ggtitle(paste0("Insert Size Distribution", "\n", infoMax, "\n", infoMin))
p1=p1 + scale_x_continuous(labels=comma)
p1=p1 + facet_grid(Library ~ Sample)
p1
dev.off()
print(warnings())

