library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
mapq=read.table(args[1], header=T)
mapq$MappingQuality=as.numeric(mapq$MappingQuality)
mapq$Count=as.numeric(mapq$Count)

# Remove leading and trailing small values
tc=sum(as.numeric(mapq$Count))
sc=sort(mapq$Count)
lowcut=sc[cumsum(sc) >= 0.01 * tc][1]
sl=unique(mapq[,c("Sample", "Library")])
mapqfilt = data.frame()
for(i in 1:nrow(sl)) {
      slc=mapq[(mapq$Sample == sl[i,"Sample"]) & (mapq$Library == sl[i,"Library"]),]
      mapqfilt=rbind(mapqfilt, slc[which(slc$Count>=lowcut)[1]:tail(which(slc$Count>=lowcut), n=1),])
}

# Plot
pdf(paste0(args[1], ".pdf"), height=4, width=8)
p1=ggplot(data=mapqfilt, aes(x=MappingQuality, y=Count))
p1=p1 + geom_line(aes(group=Library, colour=Library))
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Mapping Quality Distribution")
p1=p1 + scale_x_continuous(labels=comma, breaks=as.integer(seq(min(mapqfilt$MappingQuality), max(mapqfilt$MappingQuality), length.out=10)))
p1=p1 + facet_wrap(~ Sample) + theme(legend.position="bottom", legend.direction='vertical')
p1
print(warnings())

