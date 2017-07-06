library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
cov=read.table(args[1], header=T)
cov$Coverage=as.numeric(cov$Coverage)
cov$Count=as.numeric(cov$Count)

# Remove leading and trailing small values
tc=sum(as.numeric(cov$Count))
sc=sort(cov$Count)
lowcut=sc[cumsum(sc) >= 0.01 * tc][1]
sl=unique(cov[,c("Sample", "Library")])
covfilt = data.frame()
for(i in 1:nrow(sl)) {
      slc=cov[(cov$Sample == sl[i,"Sample"]) & (cov$Library == sl[i,"Library"]),]
      covfilt=rbind(covfilt, slc[which(slc$Count>=lowcut)[1]:tail(which(slc$Count>=lowcut), n=1),])
}
covfilt$RG = substr(covfilt$Library, 0, 12)

# Plot
pdf(paste0(args[1], ".pdf"), height=4, width=8)
p1=ggplot(data=covfilt, aes(x=Coverage, y=Count))
p1=p1 + geom_line(aes(group=RG, colour=RG))
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Coverage Distribution")
p1=p1 + scale_x_continuous(labels=comma, breaks=as.integer(seq(min(covfilt$Coverage), max(covfilt$Coverage), length.out=10)))
p1=p1 + facet_wrap(~ Sample) + theme(legend.position="bottom", legend.direction='horizontal')
p1
dev.off()
print(warnings())

