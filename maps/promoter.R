library(GenomicFeatures)

promoterTable = function(id) {
	  db = makeTxDbFromUCSC(genome=id, tablename="ccdsGene")
	  pr = reduce(promoters(db, upstream=1000, downstream=1000), ignore.strand=T)
	  pr = keepStandardChromosomes(pr)
	  pr = pr[width(pr)>1,]
	  
	  df=data.frame(chr=seqnames(pr), start=start(pr), end=end(pr), type="promoter")
	  gz = gzfile(paste0(id, ".promoter.bed.gz"), "w")
	  write.table(df, gz, quote=F, row.names=F, col.names=F, sep="\t")
	  close(gz)
}

promoterTable("hg38")
promoterTable("hg19")
promoterTable("mm10")
