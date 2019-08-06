
circJuncDEnorep<-function(files,cutoff=2,sf=T){
  
  fnames=c('chrom','start','end','name','score','strand','thickStart','thinkEnd','itemRgb','exonCount',
                       'exonSizes','exonOffsets','readNumber','circType','genename','isoformName','exonIndex/intronIndex','flankIntron')
  
  m=read.table(files[1],stringsAsFactors  = F)
  colnames(m)=fnames
  m=m[m$circType=='circRNA',]
  m=m[m$chrom!='chrM',]
  m=m[order(m$readNumber,decreasing=T),]
  g1=GRanges(seqnames = Rle(m$chrom), ranges = IRanges(m$start,m$end),strand=m$strand)
  g1$juncread=m$readNumber
  g1$gene=m$genename
  g1$exoncount=m$exonCount
  seqlevels(g1)=paste('chr',c(1:22,'X','Y'),sep='')
  

  m=read.table(files[2],stringsAsFactors  = F)
  colnames(m)=fnames
  m=m[m$circType=='circRNA',]
  m=m[m$chrom!='chrM',]
  m=m[order(m$readNumber,decreasing=T),]
  g2=GRanges(seqnames = Rle(m$chrom), ranges = IRanges(m$start,m$end),strand=m$strand)
  g2$juncread=m$readNumber
  g2$gene=m$genename
  g2$exoncount=m$exonCount
  seqlevels(g2)=paste('chr',c(1:22,'X','Y'),sep='')
  
  g1=g1[g1$juncread>=cutoff]
  g2=g2[g2$juncread>=cutoff]
  
  gg=unique(c(g1,g2))
  ngcirc=length(gg)
  
  
  juncread=matrix(0,ngcirc,2)
  colnames(juncread)=paste('juncread:',1:2,sep='')
  exoncount=matrix(0,ngcirc,2)
  colnames(exoncount)=paste('exoncount:',1:2,sep='')
  
  tmp=BiocGenerics::match(g1,gg)
  juncread[tmp,1]=g1$juncread
  exoncount[tmp,1]=g1$exoncount
  tmp=BiocGenerics::match(g2,gg)
  juncread[tmp,2]=g2$juncread
  exoncount[tmp,2]=g2$exoncount
  
  gg$juncread=juncread
  gg$exoncount=exoncount
  
  fc=rep(0,ngcirc)
  idx1=(gg$juncread[,1]!=0 & gg$juncread[,2]==0)
  idx2=(gg$juncread[,1]==0 & gg$juncread[,2]!=0)
  idx3=(gg$juncread[,1]!=0 & gg$juncread[,2]!=0)
  fc[idx1]= 0
  fc[idx2]= 10000
  if(sf){
    sfs=apply(juncread,2,sum); sfs=sfs/median(sfs)
    juncread=sweep(juncread,2,sfs,FUN='/')
  }
  fc[idx3]= juncread[idx3,2]/juncread[idx3,1]
  
  gg$fc=fc
  gg=gg[order(gg$fc,decreasing=T)]
  gg
  
}








