library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
ins=read.table(args[1], header=T)
ins$RG = substr(ins$Library, 0, 12)
upBound=max(ins[ins$Quantile >= 0.01 & ins$Quantile <= 0.99,]$InsertSize)
lowBound=min(ins[ins$Quantile >= 0.01 & ins$Quantile <= 0.99,]$InsertSize)
ins=ins[ins$InsertSize >= lowBound & ins$InsertSize <= upBound,]

# Plot
pdf(paste0(args[1], ".pdf"), height=4, width=8)
p1=ggplot(data=ins, aes(x=InsertSize, y=Count))
p1=p1 + geom_line(aes(group=Layout, colour=Layout))
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Insert Size Distribution")
p1=p1 + scale_x_continuous(labels=comma)
p1=p1 + facet_grid(Sample ~ RG)
p1=p1 + theme(axis.text.x = element_text(angle=45, hjust=1))
p1
dev.off()
print(warnings())

