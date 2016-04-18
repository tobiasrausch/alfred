library(ggplot2)
library(scales)
library(grid)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
x$status = "observed"
x$group = x$sample

# Theme
txtFontSize=10
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line.x=element_line(size=0.7, color="black"), axis.line.y=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize))

for (s in unique(x$sample)) {
    subset = x[x$sample==s,]
    totalcount = sum(as.numeric(subset$count))
    meanEst = sum(as.numeric(subset$count) * as.numeric(subset$coverage)) / totalcount
    varianceEst = sum(as.numeric(subset$coverage) * as.numeric(subset$coverage) * as.numeric(subset$count)) / totalcount - meanEst * meanEst
    sizeEst = meanEst * meanEst / (varianceEst - meanEst)
    # Poisson
    #d=data.frame(subset$coverage, dpois(subset$coverage, meanEst)*totalcount)
    # Negative Binomial
    d=data.frame(subset$coverage, dnbinom(subset$coverage, size=sizeEst, mu=meanEst) * totalcount)
    d$sample = s
    d$status = "expected"
    d$group = paste0("exp_", s)
    colnames(d)=c("coverage","count","sample", "status", "group")
    x = rbind(x, d)
    err=sum((subset$count - d$count)*(subset$count - d$count)) / totalcount
    print(s)
    print(err)
}
x$status=factor(x$status, levels=c("observed", "expected"))
x$group=factor(x$group)

# Variant allele frequency spectrum
png("coverage.png", height=800, width=1200)
p1=ggplot(data=x, aes(x=coverage, y=count)) + geom_line(aes(group=group, colour=sample, linetype=status), size=1.2) + geom_point(aes(colour=sample))
p1=p1 + xlab("Coverage Level") + ylab("#Bases at that coverage")
p1=p1 + scale_y_continuous(labels=comma) + scale_x_log10() + scienceTheme + labs(colour="Sample", shape="Status")
p1
z=dev.off()
print(warnings())

