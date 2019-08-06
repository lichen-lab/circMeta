geneExp<-function(gcirc,cufflinkfiles,bamfiles){
  
  #obtain gene coordinates
  geneexp=read.table(cufflinkfiles[1],header=T,sep="\t",stringsAsFactors = F)
  a=sapply(as.character(geneexp$locus),function(x) strsplit(x,split=":"))
  b=unlist(a)
  chrs=b[seq(length(b))%%2==1]
  locus=b[seq(length(b))%%2==0]
  d=sapply(as.character(locus),function(x) strsplit(x,split="-"))
  dd=unlist(d)
  starts=dd[seq(length(dd))%%2==1]
  ends=dd[seq(length(dd))%%2==0]
  gene.gr=GRanges(seqnames = Rle(chrs), ranges = IRanges(as.numeric(starts),as.numeric(ends)))
  gene.gr$gene=geneexp$gene_id
  ngene=length(gene.gr)
  
  fpkm=matrix(0,ngene,length(cufflinkfiles))
  for(i in 1:length(cufflinkfiles)){
    geneexp=read.table(cufflinkfiles[i],header=T,sep="\t",stringsAsFactors = F)
    id=match(gene.gr$gene,geneexp$gene_id)
    fpkm[,i]=geneexp$FPKM[id]
  }
  
  counts=getWinCounts(bamfiles,gene.gr)
  counts$fpkm=fpkm
  counts$gene.gr=gene.gr

  id=match(tolower(gcirc$gene),tolower(gene.gr$gene))
  gcirc$genecounts=counts$counts[id,]
  gcirc$genefpkm=counts$fpkm[id,]
  
  #annotate gcirc
  list(gcirc=gcirc,counts=counts)
  
}








