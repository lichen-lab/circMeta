circClass<-function(files,circ.method=c('findcirc','CIRCexplorer','CIRI'),
                    gene=NULL,gexon=NULL){
  
  circ.method=match.arg(circ.method)
  nfile=length(files)
  glist=list()
  
  #if(is.null(gene) | is.null(gexon) ){
  #  stop('gene and gexon are required!')
  #}
  
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
    fnames=circNames(circ.method)
    for(ifile in 1:nfile){
      m=read.table(files[ifile],stringsAsFactors  = F)
      colnames(m)=fnames
      m=m[m$circType=='circRNA',]
      m=m[m$chrom!='chrM',]
      m=m[order(m$readNumber,decreasing=T),]
      g=GRanges(seqnames = Rle(m$chrom), ranges = IRanges(m$start,m$end),strand=m$strand)
      g$juncread=m$readNumber
      g$gene=m$genename
      g$exoncount=m$exoncount
      glist[[ifile]]=g
      seqlevels(glist[[ifile]])=paste('chr',c(1:22,'X','Y'),sep='')
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
      g$nonjuncread=m$non_junction_reads
      g$CLR=m$junction_reads_ratio
      g$gene=m$gene_id
      glist[[ifile]]=g
      seqlevels(glist[[ifile]])=paste('chr',c(1:22,'X','Y'),sep='')
    }
  }
  
  glist=GRangesList(glist)
  gcirc=unlist(glist)
  gcirc$juncread=NULL
  gcirc.uniq=unique(gcirc)
  ngcirc=length(gcirc.uniq)
  id=BiocGenerics::match(gcirc.uniq,gcirc)
  gcirc=gcirc.uniq
  
  if(circ.method=='findcirc' | circ.method=='CIRCexplorer'){
    message('Add junction reads......')
    juncread=matrix(0,ngcirc,nfile)
    colnames(juncread)=paste('juncread:',1:nfile,sep='')
    for(ifile in 1:nfile){
      id=BiocGenerics::match(glist[[ifile]],gcirc)
      juncread[id,ifile]=glist[[ifile]]$juncread
    }
  }else if(circ.method=='CIRI'){
    message('Add junction/nonjunction reads......')
    juncread=nonjuncread=CLR=matrix(0,ngcirc,nfile)
    colnames(juncread)=paste('juncread:',1:nfile,sep='')
    colnames(nonjuncread)=paste('nonjuncread:',1:nfile,sep='')
    colnames(CLR)=paste('CLR:',1:nfile,sep='')
    for(ifile in 1:nfile){
      id=BiocGenerics::match(glist[[ifile]],gcirc)
      juncread[id,ifile]=glist[[ifile]]$juncread
      nonjuncread[id,ifile]=glist[[ifile]]$nonjuncread
      CLR[id,ifile]=glist[[ifile]]$CLR
    }
  }
  
  message('Add closest gene......')
  if(circ.method=='findcirc'){
    gcirc=addGene(gcirc,gene)
  }else if(circ.method=='CIRI'){
    gcirc$gene=gsub(',','',gcirc$gene)
  }
  
  message('Add circRNA class......')
  gcirc=addClass(gcirc,gene,gexon)
  
  if(circ.method=='findcirc' | circ.method=='CIRCexplorer'){
    sfs=apply(juncread,2,mean);sfs=sfs/median(sfs)
    juncread.norm=sweep(juncread,2,sfs,FUN='/')
    juncread.ave=apply(juncread.norm,1,mean)
    meta=data.frame(class=gcirc$class,gene=gcirc$gene,juncread,juncread.ave)
  }else if(circ.method=='CIRI'){
    CLR.ave=apply(CLR,1,mean)
    meta=data.frame(class=gcirc$class,gene=gcirc$gene,juncread,nonjuncread,CLR,CLR.ave)
  }
  
  mcols(gcirc)=meta
  gcirc
  
}





