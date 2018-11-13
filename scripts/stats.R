library(ggplot2)
library(scales)

# Output pdf
args=commandArgs(trailingOnly=TRUE)
if (length(args) > 1) { pdffile = args[2]; } else { pdffile = paste0(args[1], ".pdf"); }
pdf(pdffile)

print("Base Content")
cmd=paste0('zgrep ^BC ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=T)
for(sid in unique(all$Sample)) {
	for(rg in unique(all[all$Sample == sid,]$Library)) {
	       base = all[all$Sample == sid & all$Library == rg,]
	       p1=ggplot(data=base, aes(x=Position, y=Fraction))
	       p1=p1 + geom_line(aes(group=Base, colour=Base))
	       p1=p1 + xlab("Position in read") + ylab("Base Content Fraction")
	       p1=p1 + scale_y_continuous(labels=comma, limits=c(0, max(base$Fraction)))
	       p1=p1 + ggtitle(paste0("Base Content Distribution", "\n", "Sample: ", sid, "\n", "RG: ", rg))
	       p1=p1 + scale_x_continuous(labels=comma)
	       p1=p1 + facet_wrap(~ Read)
	       p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
               print(p1)
	}
}
print(warnings())	

print("Base Qualities")
cmd=paste0('zgrep ^BQ ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=T)
for(sid in unique(all$Sample)) {
	for(rg in unique(all[all$Sample == sid,]$Library)) {	
	       bq = all[all$Sample == sid & all$Library == rg,]
	       p1=ggplot(data=bq, aes(x=Position, y=BaseQual))
       	       p1=p1 + geom_line(aes(color=Read, group=Read))
       	       p1=p1 + xlab("Position in read") + ylab("Mean Base Quality")
       	       p1=p1 + scale_y_continuous(labels=comma, limits=c(0, max(bq$BaseQual)))
	       p1=p1 + ggtitle(paste0("Base Quality Distribution", "\n", "Sample: ", sid, "\n", "RG: ", rg))
       	       p1=p1 + scale_x_continuous(labels=comma)
       	       p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
       	       print(p1)
	}
}
print(warnings())

print("Read length")
cmd=paste0('zgrep ^RL ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=T)
for(sid in unique(all$Sample)) {
	for(rg in unique(all[all$Sample == sid,]$Library)) {
	       rl = all[all$Sample == sid & all$Library == rg,]
	       tc=sum(as.numeric(rl$Count))
	       upBound=60000
	       gr=sum(as.numeric(rl[rl$Readlength>upBound,]$Count))
	       infoMax = paste0("Read Length > ", upBound, " (", round(100 * gr / tc, digits=2), "%)")
	       rl=rl[rl$Readlength <= upBound,]
	       p1=ggplot(data=rl, aes(x=Readlength, y=Fraction))
	       p1=p1 + geom_bar(aes(fill=Read, group=Read), stat="identity", position="dodge")
	       p1=p1 + xlab("Read length") + ylab("Fraction of reads")
	       p1=p1 + scale_y_continuous(labels=comma)
	       if (gr) { p1=p1 + ggtitle(paste0("Read Length Distribution", "\n", infoMax, "\n", "Sample: ", sid, "\n", "RG: ", rg)); }
	       else { p1=p1 + ggtitle(paste0("Read Length Distribution", "\n", "Sample: ", sid, "\n", "RG: ", rg)); }
	       p1=p1 + scale_x_continuous(labels=comma)
	       p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
	       print(p1)
	}
}
print(warnings())

print("Mapping quality")
cmd=paste0('zgrep ^MQ ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=T)
for(sid in unique(all$Sample)) {
	for(rg in unique(all[all$Sample == sid,]$Library)) {
	       mq = all[all$Sample == sid & all$Library == rg,]
	       p1=ggplot(data=mq, aes(x=MappingQuality, y=Fraction))
	       p1=p1 + geom_bar(fill="darkblue", stat="identity", position="dodge")
	       p1=p1 + xlab("Mapping Quality") + ylab("Fraction of reads")
	       p1=p1 + scale_y_continuous(labels=comma)
	       p1=p1 + ggtitle(paste0("Mapping Quality Distribution", "\n", "Sample: ", sid, "\n", "RG: ", rg))
	       p1=p1 + scale_x_continuous(labels=comma)
	       p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
	       print(p1)
	}
}
print(warnings())

print("Coverage")
cmd=paste0('zgrep ^CO ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=T)
for(sid in unique(all$Sample)) {
   for(rg in unique(all[all$Sample == sid,]$Library)) {
	cov = all[all$Sample == sid & all$Library == rg,]
	tc=sum(as.numeric(cov$Count))
	if (tc > 0) {
	   upBound=max(cov[cov$Quantile >= 0.001 & cov$Quantile <= 0.999,]$Coverage)
	   gr=sum(as.numeric(cov[cov$Coverage>upBound,]$Count))
	   infoMax = paste0("Coverage > ", upBound, " (", round(100 * gr / tc, digits=2), "%)")
	   cov=cov[cov$Coverage <= upBound,]
	   p1=ggplot(data=cov, aes(x=Coverage, y=Count))
	   p1=p1 + geom_line()
	   p1=p1 + scale_y_continuous(labels=comma)
	   p1=p1 + ggtitle(paste0("Coverage Distribution", "\n", infoMax, "\n", "Sample: ", sid, "\n", "RG: ", rg))
	   p1=p1 + scale_x_continuous(labels=comma)
	   p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
	   print(p1)
	}
   }
}
print(warnings())

print("Insert Size")
cmd=paste0('zgrep ^IS ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=T)
for(sid in unique(all$Sample)) {
   for(rg in unique(all[all$Sample == sid,]$Library)) {
      ins = all[all$Sample == sid & all$Library == rg,]
      tc=sum(as.numeric(ins$Count))
      if (tc > 0) {
      	upBound=max(ins[ins$Quantile >= 0.001 & ins$Quantile <= 0.999,]$InsertSize)
	gr=sum(as.numeric(ins[ins$InsertSize>upBound,]$Count))
	infoMax = paste0("Insert size > ", upBound, " (", round(100 * gr / tc, digits=2), "%)")
	ins=ins[ins$InsertSize <= upBound,]
	p1=ggplot(data=ins, aes(x=InsertSize, y=Count))
	p1=p1 + geom_line(aes(group=Layout, colour=Layout))
	p1=p1 + scale_y_continuous(labels=comma)
	p1=p1 + ggtitle(paste0("Insert Size Distribution", "\n", infoMax, "\n", "Sample: ", sid, "\n", "RG: ", rg))
	p1=p1 + scale_x_continuous(labels=comma)
	p1=p1 + theme(axis.text.x = element_text(angle=45, hjust=1))
	print(p1)
     }
   }
}
print(warnings())

print("Homopolymer InDels")
cmd=paste0('zgrep ^IC ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=T)
for(sid in unique(all$Sample)) {
	for(rg in unique(all[all$Sample == sid,]$Library)) {	
	       ic = all[all$Sample == sid & all$Library == rg,]
	       ic$Homopolymer = factor(ic$Homopolymer, levels=c("A", "C", "G", "T", "N", "None"))
	       p1=ggplot(data=ic, aes(x=Homopolymer, y=Fraction))
       	       p1=p1 + geom_bar(aes(group=InDel, fill=InDel), stat="identity", position="dodge")
       	       p1=p1 + xlab("Homopolymer Context") + ylab("InDel Fraction")
       	       p1=p1 + scale_y_continuous(labels=comma)
	       p1=p1 + ggtitle(paste0("InDel Homopolymer Context", "\n", "Sample: ", sid, "\n", "RG: ", rg))
       	       p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
       	       print(p1)
	}
}
print(warnings())

print("InDel Size")
cmd=paste0('zgrep ^IZ ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=T)
for(sid in unique(all$Sample)) {
	for(rg in unique(all[all$Sample == sid,]$Library)) {	
	       iz = all[all$Sample == sid & all$Library == rg,]
	       p1=ggplot(data=iz, aes(x=Size, y=Count))
       	       p1=p1 + geom_line(aes(group=InDel, color=InDel))
       	       p1=p1 + xlab("InDel Size") + ylab("InDel Count")
       	       p1=p1 + scale_y_continuous(labels=comma)
      	       p1=p1 + scale_x_continuous(labels=comma)
	       p1=p1 + ggtitle(paste0("InDel Size", "\n", "Sample: ", sid, "\n", "RG: ", rg))
       	       p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
       	       print(p1)
	}
}
print(warnings())

print("GC Content")
cmd=paste0('zgrep ^GC ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=T)
if (1) {
      p1=ggplot(data=all, aes(x=GCcontent, y=fractionOfReads))
      p1=p1 + geom_line(aes(group=Library, color=Library))
      p1=p1 + xlab("GC content") + ylab("Fraction of reads or reference windows")
      p1=p1 + scale_y_continuous(labels=comma)
      p1=p1 + ggtitle("GC-Content Distribution")
      p1=p1 + scale_x_continuous(labels=comma)
      p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
      print(p1)
}
print(warnings())

# BED file
cmd=paste0('zgrep ^OT ', args[1], ' | cut -f 2-')
all=tryCatch(read.table(pipe(cmd), header=T), error=function(e) NULL)
if (!is.null(all)) {
   print("On-target rate");
   for(sid in unique(all$Sample)) {
     for(rg in unique(all[all$Sample == sid,]$Library)) {
   	ot = all[all$Sample == sid & all$Library == rg,]
	ot$Extension=as.numeric(ot$Extension)
	ot$OnTarget=as.numeric(ot$OnTarget)
	p1=ggplot(data=ot, aes(x=Extension, y=OnTarget))
	p1=p1 + geom_line()
	p1=p1 + scale_y_continuous(labels=comma, limits=c(0,1))
	p1=p1 + ggtitle(paste0("On-target rate", "\n", "Sample: ", sid, "\n", "RG: ", rg))
	p1=p1 + scale_x_continuous(labels=comma) + xlab("Left/Right Extension of target region") + ylab("Fraction of reads on-target")
	p1=p1 + facet_wrap(~ Sample) + theme(legend.position="bottom", legend.direction='horizontal')
	print(p1)
     }
   }
   print(warnings())
}
cmd=paste0('zgrep ^TC ', args[1], ' | cut -f 2-')
all=tryCatch(read.table(pipe(cmd), header=T), error=function(e) NULL)
if (!is.null(all)) {
   print("Target coverage");
   for(sid in unique(all$Sample)) {
     for(rg in unique(all[all$Sample == sid,]$Library)) {
      	x = all[all$Sample == sid & all$Library == rg,]
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
    	p1=ggplot(data=sfilt, aes(x=cov, y=frac))
    	p1=p1 + geom_line()
    	p1=p1 + scale_y_continuous(labels=comma)
	p1=p1 + ggtitle(paste0("Target coverage distribution", "\n", "Sample: ", sid, "\n", "RG: ", rg))
   	p1=p1 + scale_x_continuous(labels=comma) + xlab("Coverage") + ylab("Fraction of targets above coverage level")
    	p1=p1 + theme(legend.position="bottom", legend.direction='horizontal')
    	print(p1)
     }
   }
   print(warnings())	
}

cmd=paste0('zgrep ^PS ', args[1], ' | cut -f 2-')
all=tryCatch(read.table(pipe(cmd), header=T), error=function(e) NULL)
if (!is.null(all)) {
   print("Phased block length");
   for(sid in unique(all$Sample)) {
     for(rg in unique(all[all$Sample == sid,]$Library)) {
      	x = all[all$Sample == sid & all$Library == rg,]
    	p1=ggplot(data=x, aes(x=Size))
    	p1=p1 + geom_histogram(bins=50)
    	p1=p1 + scale_y_continuous(labels=comma)
	p1=p1 + ggtitle(paste0("Phased block length distribution", "\n", "Sample: ", sid, "\n", "RG: ", rg))
   	p1=p1 + scale_x_continuous(labels=comma) + xlab("Phased block size") + ylab("Count")
    	print(p1)
     }
   }
   print(warnings())
}


print("Metrics")
cmd=paste0('zgrep ^ME ', args[1], ' | cut -f 2-')
all=read.table(pipe(cmd), header=F, comment.char="$")
all=as.data.frame(t(all))
nc=ncol(all) + 1
me=data.frame(x=nc, y=1:nrow(all))
p1=ggplot(data=me, aes(x=x, y=y)) + geom_blank() + xlim(0, nc) + theme(line=element_blank(), text=element_blank(), title=element_blank())
for(i in 1:nrow(all)) { for(j in 1:ncol(all)) { p1=p1 + annotate("text", x=j, y=nrow(all)-i, label=all[i, j], size=2); } }
print(p1)
print(warnings())

# Close output file
dev.off()
