library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
ins=read.table(args[1], header=T)
ins$InsertSize=as.numeric(ins$InsertSize)
ins$Count=as.numeric(ins$Count)

# Remove insert size > maxvalue
tc=sum(as.numeric(ins$Count))
lastVal=max(ins$InsertSize)
gr=sum(as.numeric(ins[ins$InsertSize==lastVal,]$Count))
print(paste0("Insert size >= ", lastVal, " (", 100 * gr / tc, "%)"))
ins=ins[ins$InsertSize < lastVal,]
insfilt=ins[ins$InsertSize < 1000,]

# Remove leading and trailing small values
#sc=sort(ins$Count)
#lowcut=sc[cumsum(sc) >= 0.01 * tc][1]
#sl=unique(ins[,c("Sample", "Library")])
#insfilt = data.frame()
#for(i in 1:nrow(sl)) {
#      slc=ins[(ins$Sample == sl[i,"Sample"]) & (ins$Library == sl[i,"Library"]),]
#      insfilt=rbind(insfilt, slc[which(slc$Count>=lowcut)[1]:tail(which(slc$Count>=lowcut), n=1),])
#}

# Plot
pdf(paste0(args[1], ".pdf"))
p1=ggplot(data=insfilt, aes(x=InsertSize, y=Count))
p1=p1 + geom_line(aes(group=Layout, colour=Layout))
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Insert Size Distribution")
p1=p1 + scale_x_continuous(labels=comma, breaks=as.integer(seq(min(insfilt$InsertSize), max(insfilt$InsertSize), length.out=10)))
p1=p1 + facet_grid(Library ~ Sample)
p1
dev.off()
print(warnings())

