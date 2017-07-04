library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
base=read.table(args[1], header=T)

base$frac=rep(0,nrow(base))
sumvec=as.numeric(base[base$Base=="A",]$Count) + as.numeric(base[base$Base=="C",]$Count) + as.numeric(base[base$Base=="G",]$Count) + as.numeric(base[base$Base=="T",]$Count) + as.numeric(base[base$Base=="N",]$Count)
base[base$Base=="A",]$frac = as.numeric(base[base$Base=="A",]$Count) / as.numeric(sumvec)
base[base$Base=="C",]$frac = as.numeric(base[base$Base=="C",]$Count) / as.numeric(sumvec)
base[base$Base=="G",]$frac = as.numeric(base[base$Base=="G",]$Count) / as.numeric(sumvec)
base[base$Base=="T",]$frac = as.numeric(base[base$Base=="T",]$Count) / as.numeric(sumvec)
base[base$Base=="N",]$frac = as.numeric(base[base$Base=="N",]$Count) / as.numeric(sumvec)
d=as.data.frame(cbind(1:length(sumvec), sumvec))
colnames(d)=c("pos","Count")
tc=sum(as.numeric(d$Count))
sc=sort(d$Count)
lowcut=sc[cumsum(sc) > 0][1]
lastpos=rev(d[d$Count>=lowcut,]$pos)[1]
base=base[base$Position<lastpos,]

# Plot
png(paste0(args[1], ".png"), height=400, width=800)
p1=ggplot(data=base, aes(x=Position, y=frac))
p1=p1 + geom_line(aes(group=Base, colour=Base))
p1=p1 + xlab("Position in read") + ylab("Base Content Fraction")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Base Content Distribution")
p1=p1 + scale_x_continuous(labels=comma)
p1=p1 + facet_wrap(~ Library) + theme(legend.position="bottom", legend.direction='horizontal')
p1
dev.off()
print(warnings())

