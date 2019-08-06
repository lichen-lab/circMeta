
circFeature<-function(files,circ.method=c('findcirc','CIRCexplorer','CIRI'),
                      gene,gexon,gintron,galu,gsine,
                      atoi=NULL,DARNED=NULL,RADAR=NULL){
  
  circ.method=match.arg(circ.method)
  nfile=length(files)
  glist=list()
  
  if(is.null(gene) | is.null(gexon) | is.null(gintron) | is.null(galu) | is.null(gsine)){
    stop('gene, gexon and gintron are required!')
  }
  
  if(circ.method=='findcirc'){
    fnames=circNames(circ.method)
    for(ifile in 1:nfile){
      m=read.table(files[ifile],stringsAsFactors = F)
      colnames(m)=fnames
      m=m[m$chrom!='chrM',]
      m=m[order(m$n_reads,decreasing=T),]
      g=GRanges(seqnames = Rle(m$chrom), ranges = IRanges(m$start,m$end),strand=m$strand)
      g$juncread=m$n_reads
      glist[[ifile]]=g
      seqlevels(glist[[ifile]])=paste('chr',c(1:22,'X','Y'),sep='')
    }
  }else if(circ.method=='CIRCexplorer'){
    juncLeftIntron.list=juncRightIntron.list=list()
    fnames=circNames(circ.method)
    for(ifile in 1:nfile){
      m=read.table(files[ifile],stringsAsFactors  = F)
      colnames(m)=fnames
      m=m[m$circType=='circRNA',]
      m=m[m$chrom!='chrM',]
      m=m[order(m$readNumber,decreasing=T),]
      a=sapply(m$flankIntron,function(x) strsplit(x,split='\\|'))
      a=unlist(a)
      leftIntron=a[seq(length(a))%%2==1]
      rightIntron=a[seq(length(a))%%2==0]
      id1=leftIntron!="None" & leftIntron!=''
      id2=rightIntron!="None" & rightIntron!=''
      id=id1 & id2
      g=GRanges(seqnames = Rle(m$chrom), ranges = IRanges(m$start,m$end),strand=m$strand)
      g$juncread=m$readNumber
      g$gene=m$genename
      g$exoncount=m$exoncount
      glist[[ifile]]=g[which(id)]
      seqlevels(glist[[ifile]])=paste('chr',c(1:22,'X','Y'),sep='')
      juncLeftIntron.list[[ifile]]=stringToGRange(leftIntron[which(id)])
      juncRightIntron.list[[ifile]]=stringToGRange(rightIntron[which(id)])
    }
  }else if(circ.method=='CIRI'){
    fnames=circNames(circ.method)
    for(ifile in 1:nfile){
      m=read.table(files[ifile],stringsAsFactors  = F,sep='\t',skip=1)
      m=m[,-12]
      colnames(m)=fnames[-12]
      m=m[m$chr!='chrM',]
      m=m[order(m$junction_reads,decreasing=T),]
      g=GRanges(seqnames = Rle(m$chr), ranges = IRanges(m$circRNA_start,m$circRNA_end),strand=m$strand)
      g$juncread=m$junction_reads
      g$gene=m$gene_id
      g$juncread=m$junction_reads
      g$nonjuncread=m$non_junction_reads
      g$CLR=m$junction_reads_ratio
      glist[[ifile]]=g
      seqlevels(glist[[ifile]])=paste('chr',c(1:22,'X','Y'),sep='')
    }
  }
  
  glist=GRangesList(glist)
  gcirc=unlist(glist)
  gcirc$juncread=NULL
  gcirc.uniq=unique(gcirc)
  ngcirc=length(gcirc.uniq)
  id.uniq=BiocGenerics::match(gcirc.uniq,gcirc)
  gcirc=gcirc.uniq
  
  if(circ.method=='findcirc' | circ.method=='CIRCexplorer'){
    message('Add junction reads......')
    juncread=matrix(0,ngcirc,nfile)
    colnames(juncread)=paste('juncread:',1:nfile,sep='')
    for(ifile in 1:nfile){
      tmp=BiocGenerics::match(glist[[ifile]],gcirc)
      juncread[tmp,ifile]=glist[[ifile]]$juncread
    }
  }else if(circ.method=='CIRI'){
    message('Add junction/nonjunction reads......')
    juncread=nonjuncread=CLR=matrix(0,ngcirc,nfile)
    colnames(juncread)=paste('juncread:',1:nfile,sep='')
    colnames(nonjuncread)=paste('nonjuncread:',1:nfile,sep='')
    colnames(CLR)=paste('CLR:',1:nfile,sep='')
    for(ifile in 1:nfile){
      tmp=BiocGenerics::match(glist[[ifile]],gcirc)
      juncread[tmp,ifile]=glist[[ifile]]$juncread
      nonjuncread[tmp,ifile]=glist[[ifile]]$nonjuncread
      CLR[tmp,ifile]=glist[[ifile]]$CLR
    }
  }
  
  if(circ.method=='findcirc' | circ.method=='CIRI'){
    res=getJuncIntron(gcirc,gintron,gexon)
    tmp=BiocGenerics::match(res$gcirc,gcirc)
    juncread=juncread[tmp,]
    if(circ.method=='CIRI'){
      nonjuncread=nonjuncread[tmp,]
      CLR=CLR[tmp,]
    }
    gcirc=res$gcirc
    juncLeftIntron=res$juncLeftIntron
    juncRightIntron=res$juncRightIntron
  }else if(circ.method=='CIRCexplorer'){
    juncLeftIntron.list=GRangesList(juncLeftIntron.list)
    juncRightIntron.list=GRangesList(juncRightIntron.list)
    juncLeftIntron=unlist(juncLeftIntron.list)
    juncRightIntron=unlist(juncRightIntron.list)
    juncLeftIntron=juncLeftIntron[id.uniq]
    juncRightIntron=juncRightIntron[id.uniq]
  }
  
  message('Add closest gene......')
  if(circ.method=='findcirc'){
    gcirc=addGene(gcirc,gene)
  }else if(circ.method=='CIRI'){
    gcirc$gene=gsub(',','',gcirc$gene)
  }
  
  message('Add circRNA class......')
  gcirc=addClass(gcirc,gene,gexon)
  
  message('Add intron length......')
  intronlen=intronLen(juncLeftIntron,juncRightIntron)
  
  message('Add Alu enrichment......')
  aluenrich=aluEnrich(galu,juncLeftIntron,juncRightIntron)
  message('Add Alu score......')
  aluscore=aluScore(gcirc,galu,juncLeftIntron,juncRightIntron)
  message('Add SINE enrichment......')
  SINEenrich=SINEEnrich(gsine,juncLeftIntron,juncRightIntron)
  message('Add SINE score......')
  SINEscore=SINEScore(gcirc,gsine,juncLeftIntron,juncRightIntron)
  
  if(!is.null(atoi)){
    message('Add AtoI enrichment in intron......')
    atoienrichintron=atoiEnrichIntron(atoi,juncLeftIntron,juncRightIntron)
    message('Add AtoI enrichment in Alu')
    atoienrichalu=atoiEnrichAlu(gcirc,atoi,galu,juncLeftIntron,juncRightIntron)
  }
  
  if(circ.method=='findcirc' | circ.method=='CIRCexplorer'){
    sfs=apply(juncread,2,mean);sfs=sfs/median(sfs)
    juncread.norm=sweep(juncread,2,sfs,FUN='/')
    juncread.ave=apply(juncread.norm,1,mean)
    meta=data.frame(class=gcirc$class,gene=gcirc$gene,intronlen,aluenrich,aluscore,SINEenrich,SINEscore,juncread,juncread.ave)
  }else if(circ.method=='CIRI'){
    CLR.ave=apply(CLR,1,mean)
    meta=data.frame(class=gcirc$class,gene=gcirc$gene,intronlen,aluenrich,aluscore,SINEenrich,SINEscore,juncread,nonjuncread,CLR,CLR.ave)
  }
  
  if(!is.null(atoi)){
    meta$atoienrichintron=atoienrichintron
    meta$atoienrichalu=atoienrichalu
  }
  mcols(gcirc)=meta
  gcirc$juncLeftIntron.start=start(juncLeftIntron)
  gcirc$juncLeftIntron.end=end(juncLeftIntron)
  gcirc$juncRightIntron.start=start(juncRightIntron)
  gcirc$juncRightIntron.end=end(juncRightIntron)
  
  gcirc
  #list(juncRightIntron=juncRightIntron,juncLeftIntron=juncLeftIntron,gcirc=gcirc)
  
}