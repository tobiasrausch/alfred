library(ggplot2)
library(scales)

args=commandArgs(trailingOnly=TRUE)
x=read.table(args[1], header=T)

s = data.frame()
sl=unique(x[,c("Sample", "Library")])
for(i in 1:nrow(sl)) {
      slc=x[(x$Sample == sl[i,"Sample"]) & (x$Library == sl[i,"Library"]),]
      for (j in 0:max(x$AvgCov)) { s=rbind(s, cbind(Sample=as.character(sl[i,"Sample"]), Library=as.character(sl[i,"Library"]), frac=mean(x$AvgCov >= j), cov=j)); }
}
s$cov = as.numeric(as.character(s$cov))
s$frac = as.numeric(as.character(s$frac))

# Remove trailing small values
tc=sum(as.numeric(s$frac))
sc=sort(s$frac)
lowcut=sc[cumsum(sc) >= 0.01 * tc][1]
sl=unique(s[,c("Sample", "Library")])
sfilt = data.frame()
for(i in 1:nrow(sl)) {
      slc=s[(s$Sample == sl[i,"Sample"]) & (s$Library == sl[i,"Library"]),]
      sfilt=rbind(sfilt, slc[which(slc$frac>=lowcut)[1]:tail(which(slc$frac>=lowcut), n=1),])
}

# Plot
pdf(paste0(args[1], ".pdf"), height=4, width=8)
p1=ggplot(data=sfilt, aes(x=cov, y=frac))
p1=p1 + geom_line(aes(group=Library, colour=Library))
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Target coverage distribution")
p1=p1 + scale_x_continuous(labels=comma) + xlab("Coverage") + ylab("Fraction of targets above coverage level")
p1=p1 + facet_wrap( ~ Sample) + theme(legend.position="bottom", legend.direction='horizontal')
p1
dev.off()
print(warnings())

