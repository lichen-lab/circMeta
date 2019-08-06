circJuncDE<-function(files,designs,
                 circ.method=c('findcirc','CIRCexplorer','CIRI'),
                 DE.method=c('pois.ztest','DESeq2','edgeR','pois.glm','nb.glm'),
                 gene=NULL,gexon=NULL,
                 cutoff=2,sf=TRUE){
  
  circ.method=match.arg(circ.method)
  DE.method=match.arg(DE.method)
  nfile=length(files)
  glist=list()
  
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
  #id.uniq=match(gcirc.uniq,gcirc)
  gcirc=gcirc.uniq
  
  message('Add closest gene......')
  if(circ.method=='findcirc'){
    if(!is.null(gene)) gcirc=addGene(gcirc,gene)
  }else if(circ.method=='CIRI'){
    gcirc$gene=gsub(',','',gcirc$gene)
  }
  
  message('Add circRNA class......')
  if(!is.null(gene) & !is.null(gexon)) gcirc=addClass(gcirc,gene,gexon)
  
  if(circ.method=='findcirc' | circ.method=='CIRCexplorer'){
    message('Add junction reads......')
    juncread=matrix(0,ngcirc,nfile)
    colnames(juncread)=paste('juncread:',designs,1:nfile,sep='')
    for(ifile in 1:nfile){
      tmp=BiocGenerics::match(glist[[ifile]],gcirc)
      juncread[tmp,ifile]=glist[[ifile]]$juncread
    }
  }else if(circ.method=='CIRI'){
    message('Add junction/nonjunction reads......')
    juncread=nonjuncread=CLR=matrix(0,ngcirc,nfile)
    colnames(juncread)=paste('juncread:',designs,1:nfile,sep='')
    colnames(nonjuncread)=paste('nonjuncread:',designs,1:nfile,sep='')
    colnames(CLR)=paste('CLR:',designs,1:nfile,sep='')
    for(ifile in 1:nfile){
      tmp=BiocGenerics::match(glist[[ifile]],gcirc)
      juncread[tmp,ifile]=glist[[ifile]]$juncread
      nonjuncread[tmp,ifile]=glist[[ifile]]$nonjuncread
      CLR[tmp,ifile]=glist[[ifile]]$CLR
    }
  }
  
  
  dat=NULL
  dat$counts=juncread
  dat$gcirc=gcirc
  dat$designs=designs
  dat$cutoff=cutoff
  dat$sf=sf
  dat$circ.method=circ.method
  if(circ.method=='CIRI'){
    dat$nonjuncread=nonjuncread
    dat$CLR=CLR
  }
  if(DE.method=='DESeq2') res=runDESeq2(dat)
  if(DE.method=='edgeR') res=runedgeR(dat)
  if(DE.method=='pois.glm') res=runPois.glm(dat)
  if(DE.method=='pois.ztest') res=runPois.ztest(dat)
  if(DE.method=='nb.glm') res=runNB.glm(dat)
  
  res
  
}





