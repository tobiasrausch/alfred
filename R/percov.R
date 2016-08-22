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

print(summary(s))
# Plot
pdf(paste0(args[1], ".pdf"))
p1=ggplot(data=s, aes(x=cov, y=frac))
p1=p1 + geom_line(aes(group=Library, colour=Library))
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("Target coverage distribution")
p1=p1 + scale_x_continuous(labels=comma) + xlab("Coverage") + ylab("Fraction of targets above coverage level")
p1=p1 + facet_wrap( ~ Sample) + theme(legend.position="bottom", legend.direction='vertical')
p1
print(warnings())

