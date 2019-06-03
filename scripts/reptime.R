library(reshape2)
library(ggplot2)
library(scales)
library(optparse)

optl = list(
     make_option(c("-f", "--file"), type="character", default=NULL, help="input profile file"),
     make_option(c("-r", "--region"), type="character", default=NULL, help="region chr:start-end")
     )

opt_parser = OptionParser(option_list=optl)
opt = parse_args(opt_parser)

if ((is.null(opt$file)) || (is.null(opt$region))) {
   print_help(opt_parser)
   stop("Input options missing!", call.=FALSE)
}

x = read.table(opt$file, header=T)

if (!is.null(opt$region)) {
   s=unlist(strsplit(opt$region, ":"))
   chr=s[1]
   minStart = 0
   maxStart = max(x[x$chr==chr,]$pos) + 1
   if (length(s) > 1) {
      v=unlist(strsplit(s[2], "-"))
      minStart = as.integer(v[1])
      maxStart = as.integer(v[2])
   }
}

x$grp = cumsum(x$reptime == -1)
x = x[x$reptime >= 0,]
x$grp = factor(x$grp)
x = x[x$chr == chr & x$pos>=minStart & x$pos<=maxStart,];
x$chr = factor(x$chr)


p1 = ggplot(data=x, aes(x=pos, y=reptime))
p1 = p1 + geom_line(aes(group=grp))
p1 = p1 + xlab(paste0(chr, " position")) + ylab("Replication Time")
p1 = p1 + scale_x_continuous(labels=comma)
p1 = p1 + scale_y_continuous(labels=comma)
ggsave(paste0(chr, ".", minStart, ".", maxStart ,".reptime.png"), width=16, height=4)
print(warnings())

