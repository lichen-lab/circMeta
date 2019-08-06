
####################################
### circFeature/circClass module

stringToGRange<-function(string){
  a=sapply(as.character(string),function(x) strsplit(x,split=":"))
  b=unlist(a)
  chrs=b[seq(length(b))%%2==1]
  locus=b[seq(length(b))%%2==0]
  d=sapply(as.character(locus),function(x) strsplit(x,split="-"))
  dd=unlist(d)
  starts=dd[seq(length(dd))%%2==1]
  ends=dd[seq(length(dd))%%2==0]
  starts=as.numeric(starts)
  ends=as.numeric(ends)
  g=GRanges(seqnames = Rle(chrs), ranges = IRanges(starts,ends))
  g
}


circNames<-function(circ.method=c('findcirc','CIRCexplorer','CIRI')){
  circ.method=match.arg(circ.method)
  findcirc.names=c('chrom','start','end','name','n_reads','strand','n_uniq','uniqbridges','best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap','breakpoints','signal','strandmatch','category')
  
  circexplorer.names=c('chrom','start','end','name','score','strand','thickStart','thinkEnd','itemRgb','exonCount',
                       'exonSizes','exonOffsets','readNumber','circType','genename','isoformName','exonIndex/intronIndex','flankIntron')
  
  ciri.names=c('circRNA_ID','chr','circRNA_start','circRNA_end','junction_reads','SM_MS_SMS','non_junction_reads','junction_reads_ratio','circRNA_type','gene_id','strand','junction_reads_ID')
  
  if(circ.method=='findcirc'){
    return(findcirc.names)
  }else if(circ.method=='CIRCexplorer'){
    return(circexplorer.names)
  }else if(circ.method=='CIRI'){
    return(ciri.names)
  }
}




addGene<-function(gcirc,gene){
  id=nearest(gcirc,gene)
  gcirc$gene=gene$genename[id]
  gcirc
}


addClass<-function(gcirc,gene,gexon){
  ngcirc=length(gcirc)
  gclass=character(ngcirc)
  g1=g2=gcirc
  end(g1)=start(g1)
  start(g2)=end(g2)
  o1=countOverlaps(g1,gene)
  o2=countOverlaps(g2,gene)
  id1=which(o1==0 & o2==0)
  gclass[id1]='intergenetic'
  id2=setdiff(seq(ngcirc),id1)
  o3=countOverlaps(gcirc[id2],gexon)
  gclass[id2[o3>0]]='exon'
  gclass[id2[o3==0]]='intron'
  gcirc$class=gclass
  gcirc
}


getJuncIntron<-function(gcirc,gintron,gexon){
  ncirc=length(gcirc)
  nintron=length(gintron)
  nexon=length(gexon)
  
  ## obtain circ with at least one exon
  gexon=sort(gexon)
  o=findOverlaps(gcirc,gexon)
  oo=as.data.frame(o)
  colnames(oo)=c('q','s')
  id.noexon=setdiff(seq(ncirc),unique(oo$q))
  id.exon=setdiff(seq(ncirc),id.noexon)
  gcirc.exon=gcirc[id.exon]
  ngcirc.exon=length(gcirc.exon)
  
  ## for circ with at least one exon, exon on leftest and rightest
  o=findOverlaps(gcirc.exon,gexon)
  oo=as.data.frame(o)
  exon.min=tapply(oo$subjectHits,oo$queryHits, min)
  exon.min=gexon[exon.min]
  exon.max=tapply(oo$subjectHits,oo$queryHits, max)
  exon.max=gexon[exon.max]
  
  #splicing site left, find exon-intron junctions
  exon.min.left=exon.min
  intron.left=exon.min
  end(exon.min.left)=start(exon.min.left)
  o1=findOverlaps(exon.min.left,gintron)
  oo1=as.data.frame(o1)
  colnames(oo1)=c('q','s')
  dis=end(gintron[oo1$s])-start(exon.min.left[oo1$q])
  oo1=oo1[(abs(dis)==0 | abs(dis)==1),]
  intron.tmp=gintron[oo1$s]
  len.tmp=end(intron.tmp)-start(intron.tmp)
  tt=tapply(len.tmp,oo1$q, function(x){which.min(x)})
  tb=table(oo1$q)
  tbs=c(cumsum(tb))
  left=sapply(1:length(tt),function(i) tbs[i-1]+tt[i])
  left=unlist(left)
  left=c(1,left)
  intron.left=intron.tmp[left]
  
  #splicing site right, find exon-intron junctions
  exon.max.right=exon.max
  intron.right=exon.max
  start(exon.max.right)=end(exon.max.right)
  o2=findOverlaps(exon.max.right,gintron)
  oo2=as.data.frame(o2)
  colnames(oo2)=c('q','s')
  dis=start(gintron[oo2$s])-end(exon.max.right[oo2$q])
  oo2=oo2[(abs(dis)==0 | abs(dis)==1),]
  intron.tmp=gintron[oo2$s]
  len.tmp=end(intron.tmp)-start(intron.tmp)
  tt=tapply(len.tmp,oo2$q, function(x){which.min(x)})
  tb=table(oo2$q)
  tbs=c(cumsum(tb))
  right=sapply(1:length(tt),function(i) tbs[i-1]+tt[i])
  right=unlist(right)
  right=c(1,right)
  intron.right=intron.tmp[right]
  
  #pair
  id=intersect(unique(oo1$q),unique(oo2$q))
  gcirc.exon=gcirc.exon[id]
  intron.left=intron.left[match(id,unique(oo1$q))]
  intron.right=intron.right[match(id,unique(oo2$q))]
  d1=(end(intron.left)-start(gcirc.exon))
  d2=(start(intron.right)-end(gcirc.exon))
  id=(abs(d1)==0 | abs(d1)==1) & (abs(d2)==1 | abs(d2)==0)
  juncRightIntron=intron.right[id]
  juncLeftIntron=intron.left[id]
  gcirc=gcirc.exon[id]
  exon.juncLeftIntron=exon.min
  exon.juncRightIntron=exon.max
  
  list(juncRightIntron=juncRightIntron,juncLeftIntron=juncLeftIntron,
       gcirc=gcirc)
}



intronLen<-function(juncLeftIntron,juncRightIntron){
  lenLeft=end(juncLeftIntron)-start(juncLeftIntron)
  lenRight=end(juncRightIntron)-start(juncRightIntron)
  len=round((lenLeft+lenRight)/2)
  len
}


aluEnrich<-function(galu,juncLeftIntron,juncRightIntron){
  o1=findOverlaps(juncLeftIntron,galu)
  tb1=table(queryHits(o1))
  o2=findOverlaps(juncRightIntron,galu)
  tb2=table(queryHits(o2))
  tb.names=unique(c(names(tb1),names(tb2)))
  tb=numeric(length(tb.names)); names(tb)=tb.names
  tb[match(names(tb1),tb.names)]=tb1
  tb[match(names(tb2),tb.names)]=tb2+tb[match(names(tb2),tb.names)]
  a=numeric(length(juncLeftIntron))
  a[as.numeric(names(tb))]=tb
  a
}



aluScore<-function(gcirc,galu,juncLeftIntron,juncRightIntron){
  o1=findOverlaps(juncLeftIntron,galu,ignore.strand=T)
  id1=unique(queryHits(o1))
  o2=findOverlaps(juncRightIntron,galu,ignore.strand=T)
  id2=unique(queryHits(o2))
  ids=intersect(id1,id2)
  o1=as.data.frame(o1)
  o2=as.data.frame(o2)
  colnames(o1)=colnames(o2)=c('q','s')
  o1=o1[o1$q %in% ids,]
  o2=o2[o2$q %in% ids,]
  strand=as.data.frame(galu)[,5]
  o1=data.frame(q=o1$q,s=o1$s,strand=strand[o1$s])
  #tb1=tapply(o1$strand,o1$q, function(x){ tb=table(x);tb[1]*tb[2]})
  o2=data.frame(q=o2$q,s=o2$s,strand=strand[o2$s])
  #tb2=tapply(o2$strand,o2$q, function(x){ tb=table(x);tb[1]*tb[2]})
  tb1=tapply(o1$strand,o1$q, function(x){ table(x)[1:2]})
  tb1=matrix(unlist(tb1),length(tb1),2,byrow=T)
  tb2=tapply(o2$strand,o2$q, function(x){ table(x)[1:2]})
  tb2=matrix(unlist(tb2),length(tb2),2,byrow=T)
  #IRAlusacross – IRAluswithin
  score1=tb1[,1]*tb1[,2]
  score2=tb2[,1]*tb2[,2]
  score3=tb1[,1]*tb2[,2]+tb1[,2]*tb2[,1]
  score=score3-apply(cbind(score1,score2),1,max)
  score[score<0]=0
  aluscore=numeric(length(gcirc))
  aluscore[ids]=score
  aluscore
  
}



SINEEnrich<-function(gsine,juncLeftIntron,juncRightIntron){
  o1=findOverlaps(juncLeftIntron,gsine)
  tb1=table(queryHits(o1))
  o2=findOverlaps(juncRightIntron,gsine)
  tb2=table(queryHits(o2))
  tb.names=unique(c(names(tb1),names(tb2)))
  tb=numeric(length(tb.names)); names(tb)=tb.names
  tb[match(names(tb1),tb.names)]=tb1
  tb[match(names(tb2),tb.names)]=tb2+tb[match(names(tb2),tb.names)]
  a=numeric(length(juncLeftIntron))
  a[as.numeric(names(tb))]=tb
  a
}



SINEScore<-function(gcirc,gsine,juncLeftIntron,juncRightIntron){
  o1=findOverlaps(juncLeftIntron,gsine,ignore.strand=T)
  id1=unique(queryHits(o1))
  o2=findOverlaps(juncRightIntron,gsine,ignore.strand=T)
  id2=unique(queryHits(o2))
  ids=intersect(id1,id2)
  o1=as.data.frame(o1)
  o2=as.data.frame(o2)
  colnames(o1)=colnames(o2)=c('q','s')
  o1=o1[o1$q %in% ids,]
  o2=o2[o2$q %in% ids,]
  strand=as.data.frame(gsine)[,5]
  o1=data.frame(q=o1$q,s=o1$s,strand=strand[o1$s])
  #tb1=tapply(o1$strand,o1$q, function(x){ tb=table(x);tb[1]*tb[2]})
  o2=data.frame(q=o2$q,s=o2$s,strand=strand[o2$s])
  #tb2=tapply(o2$strand,o2$q, function(x){ tb=table(x);tb[1]*tb[2]})
  tb1=tapply(o1$strand,o1$q, function(x){ table(x)[1:2]})
  tb1=matrix(unlist(tb1),length(tb1),2,byrow=T)
  tb2=tapply(o2$strand,o2$q, function(x){ table(x)[1:2]})
  tb2=matrix(unlist(tb2),length(tb2),2,byrow=T)
  #IRSINEsacross – IRSINEswithin
  score1=tb1[,1]*tb1[,2]
  score2=tb2[,1]*tb2[,2]
  score3=tb1[,1]*tb2[,2]+tb1[,2]*tb2[,1]
  score=score3-apply(cbind(score1,score2),1,max)
  score[score<0]=0
  SINEscore=numeric(length(gcirc))
  SINEscore[ids]=score
  SINEscore
}


getAtoI<-function(vcf,DARNED=NULL,RADAR=NULL){
  atoi=read.table(vcf,stringsAsFactors=F)
  atoi=atoi[atoi[,3]=='A' & atoi[,4]=='G',]
  atoi=GRanges (seqnames = Rle(atoi[,1]), ranges = IRanges(atoi[,2],atoi[,2]))
  if(!is.null(DARNED)){
    o1=countOverlaps(atoi,DARNED)
    atoi=atoi[o1>0]
  }else if(!is.null(RADAR)){
    o2=countOverlaps(atoi,RADAR)
    atoi=atoi[o2>0]
  }else if(!is.null(DARNED) & !is.null(RADAR)){
    o1=countOverlaps(atoi,DARNED)
    o2=countOverlaps(atoi,RADAR)
    atoi=atoi[o1>0 & o2>0]
  }
  atoi
}


atoiEnrichIntron<-function(atoi,juncLeftIntron,juncRightIntron){
  atoi.left=countOverlaps(juncLeftIntron,atoi)
  atoi.right=countOverlaps(juncRightIntron,atoi)
  (atoi.left+atoi.right)/2
}


atoiEnrichAlu<-function(gcirc,atoi,galu,juncLeftIntron,juncRightIntron){
  o1=findOverlaps(juncLeftIntron,galu)
  oo1=as.data.frame(o1)
  colnames(oo1)=c('q','s')
  c1=countOverlaps(galu[oo1$s],atoi)
  t1=tapply(c1,oo1$q,sum) 
  num1=numeric(length(gcirc))
  num1[unique(oo1$q)]=t1
  o2=findOverlaps(juncRightIntron,galu)
  oo2=as.data.frame(o2)
  colnames(oo2)=c('q','s')
  c2=countOverlaps(galu[oo2$s],atoi)
  t2=tapply(c2,oo2$q,sum) 
  num2=numeric(length(gcirc))
  num2[unique(oo2$q)]=t2
  num=(num1+num2)/2
  num
}




####################################
### circJoint module

read.BAM<-function (fn, ext = 0){
  what = c("rname", "strand", "pos", "qwidth")
  TSS.counts = NULL
  param = ScanBamParam(what = what)
  bam = scanBam(fn, param = param)[[1]]
  ix = !is.na(bam$rname) & !is.na(bam$pos)
  if (ext == 0) {
    qwidth = bam$qwidth[ix]
    IRange.reads = GRanges(seqnames = Rle(bam$rname[ix]),
                           ranges = IRanges(bam$pos[ix], width = bam$qwidth[ix]))
  }
  else {
    idx.minus = bam$strand == "-"
    idx.minus = idx.minus[ix]
    ss = bam$pos[ix]
    ee = bam$pos[ix] + bam$qwidth[ix]
    ss2 = ss
    ss2[idx.minus] = ee[idx.minus] - ext
    IRange.reads = GRanges(seqnames = Rle(bam$rname[ix]),
                           ranges = IRanges(ss2, width = ext))
  }
  IRange.reads
}



getWinCounts<-function (bamfiles, wins)
{
  if (class(wins) == "data.frame") {
    ran = IRanges(start = wins$start, end = wins$end)
    wins.ranges = GRanges(seqnames = Rle(wins$chr), ranges = ran)
  }
  else if (class(wins) == "GRanges") {
    wins.ranges = wins
  }
  else stop("input genomic window must be a GRanges or a data frame!")
  result = matrix(0, nrow = length(wins.ranges), ncol = length(bamfiles))
  totalCounts = rep(0, length(bamfiles))
  for (i in 1:length(bamfiles)) {
    reads.ranges = read.BAM(bamfiles[i])
    result[, i] = countOverlaps(wins.ranges, reads.ranges)
    totalCounts[i] = length(reads.ranges)
  }
  names(totalCounts) = colnames(result) = bamfiles
  list(counts = result, totalCounts = totalCounts)
}



####################################
### circCLRDE

rowVars<-function (x, center=NULL, ...) {
  n=!is.na(x)
  n=rowSums(n)
  n[n <= 1]=NA
  if (is.null(center)) {
    center=rowMeans(x, ...)
  }
  x=x - center
  x=x * x
  x=rowSums(x, ...)
  x=x/(n - 1)
  x
}

estRho<-function(X, N) {
  p=X/N
  mu=rowMeans(p)
  mu[mu==0]=1e-5
  mu[mu==1]=1-1e-5
  vv=rowVars(p)
  rho=vv/mu/(1-mu)
  rho[rho>=1]=1-1e-5
  rho[rho==0]=1e-5
  rho
}


estMu<-function(X, N) {
  p=X/N
  const=mean(p, na.rm=TRUE)
  mu=(rowSums(X)+const)/(rowSums(N)+1)
  nrep=ncol(X)
  matrix(rep(mu,nrep), ncol=nrep)
}



estRhoBayes<-function(X, N, est.mu) {
  prior=estRhoPrior(X, N)
  estRhoShrink(X, N, prior, est.mu)
}


estRhoPrior<-function(X, N) {
  if(ncol(X) == 1) # single rep
    return(c(-3, 1))
  
  ## keep circRNA with large read counts
  ix=rowMeans(N>10)==1 & rowSums(N==0)==0
  if(sum(ix) < 50) {
    warning("The coverages are too low. Cannot get good estimations of prior. Use arbitrary prior N(-3,1).")
    return(c(-3, 1))
  }
  
  X=X[ix,,drop=FALSE]; N=N[ix,,drop=FALSE]
  p=X/N
  mu=rowMeans(p)
  mu[mu==0]=1e-5
  mu[mu==1]=1-1e-5
  vv=rowVars(p)
  rho=vv/mu/(1-mu)
  rho=rho[vv>0]
  lrho=log(rho[rho>0])
  prior.mean=median(lrho, na.rm=TRUE)
  prior.sd=IQR(lrho, na.rm=TRUE) /1.349
  c(prior.mean, prior.sd)
}



estRhoShrink<-function(X, N, prior, est.mu) {
  ## penalized likelihood function
  plik.logN<-function(n,x,mu,m0,tau,rho)
    -(sum(dbb(n,x, mu,exp(rho)))+dnorm(rho, mean=m0, sd=tau, log=TRUE))
  
  if(!is.matrix(est.mu))
    est.mu=as.matrix(est.mu)
  
  shrk.rho=exp(rep(prior[1],nrow(N)))
  ix=rowSums(N>0) > 0
  X2=X[ix, ,drop=FALSE]; N2=N[ix,,drop=FALSE]; est.mu2=est.mu[ix,,drop=FALSE]
  shrk.rho2=rep(0, nrow(X2))
  for(i in 1:nrow(X2)) {
    shrk.one=optimize(f=plik.logN, n=N2[i,], x=X2[i,], mu=est.mu2[i,], m0=prior[1], tau=prior[2],
                      interval=c(-5, log(0.99)),tol=1e-3)
    shrk.rho2[i]=exp(shrk.one$minimum)
  }
  shrk.rho[ix]=shrk.rho2
  shrk.rho
}


dbb<-function(n, x, mu, rho, log=TRUE)  {
  tmp=1/rho-1
  alpha=mu*tmp
  beta=tmp - alpha
  #lbeta B(a,b)=Γ(a)Γ(b)/Γ(a+b) on log
  #lchoose choose(n, k). For k ≥ 1 it is defined as n(n-1)…(n-k+1) / k!
  #https://en.wikipedia.org/wiki/Beta-binomial_distribution
  v=lchoose(n,x)-lbeta(beta, alpha)+lbeta(n-x + beta,x+alpha)
  if(!log)
    return(exp(v))
  else return(v)
}


rbb<-function(n, mu, rho) {
  tmp=1/rho-1
  alpha=mu*tmp
  beta=tmp - alpha
  a=rbeta(length(mu),alpha,beta)
  x=rbinom(length(mu),n,a)
  x
}


waldTest<-function(n1, est.mu1, rho1, n2, est.mu2, rho2) {
  dif=est.mu1 - est.mu2
  n1m=rowSums(n1)  
  n2m=rowSums(n2)
  var1=rowSums(n1*est.mu1*(1-est.mu1)*(1+(n1-1)*rho1)) / (n1m)^2
  var2=rowSums(n2*est.mu2*(1-est.mu2)*(1+(n2-1)*rho2)) / (n2m)^2
  vv=var1 + var2
  vv[vv<1e-5]=1e-5
  se=sqrt(vv)
  stat=dif/se
  pval=2 * pnorm(-abs(stat))
  fdr=p.adjust(pval, method="fdr")
  data.frame(stat=stat,pval=pval, fdr=fdr)
}


lrTest<-function(x1, n1, est.mu1, rho1, x2, n2, est.mu2, rho2,est.mu) {
  l0=rowSums(dbb(n1,x1,est.mu,rho1)+dbb(n2,x2,est.mu,rho2) )
  l1=rowSums(dbb(n1,x1,est.mu1,rho1)+dbb(n2,x2,est.mu2,rho2) )
  lr=2*(l1-l0)
  pval=1-pchisq(lr, df=1)
  fdr=p.adjust(pval, method="fdr")
  data.frame(stat=lr,pval=pval, fdr=fdr)
}


fisherTest<-function(x1,n1,x2,n2){
  xx=cbind(rowSums(x1),rowSums(n1),rowSums(x2),rowSums(n2))
  pval=apply(xx,1,function(xxx) {fisher.test(matrix(c(xxx[1],xxx[2]-xxx[1],xxx[3],xxx[4]-xxx[3]),2,2,byrow=T))$p.value})
  fdr=p.adjust(pval, method="fdr")
  data.frame(pval=pval, fdr=fdr)
}


chisqTest<-function(x1,n1,x2,n2){
  xx=cbind(rowSums(x1),rowSums(n1),rowSums(x2),rowSums(n2))
  pval=apply(xx,1,function(xxx) {chisq.test(matrix(c(xxx[1],xxx[2]-xxx[1]+1,xxx[3],xxx[4]-xxx[3]+1),2,2,byrow=T))$p.value})
  fdr=p.adjust(pval, method="fdr")
  data.frame(pval=pval, fdr=fdr)
}



getCLR<-function(res){
  x1=res[,grep('^juncread:0',colnames(res))]
  n1=res[,grep('^nonjuncread:0',colnames(res))]
  n1=n1+x1
  x2=res[,grep('^juncread:1',colnames(res))]
  n2=res[,grep('^nonjuncread:1',colnames(res))]
  n2=n2+x2
  x1=as.matrix(x1);x2=as.matrix(x2)
  n1=as.matrix(n1);n2=as.matrix(n2)
  list(x1=x1,n1=n1,x2=x2,n2=n2)
}


#######################################
### DE for circRNA (used in simulation)

.runedgeR <- function(dat) {
  require(edgeR)
  d=DGEList(counts=dat$counts, group=dat$designs)
  if(dat$sf){
    d=calcNormFactors(d)
  }else{
    d=calcNormFactors(d,method='none')
  }
  d=estimateCommonDisp(d)
  d=estimateTagwiseDisp(d)
  fit=exactTest(d)
  res=topTags(fit, n=nrow(dat$counts))
  res=as.data.frame(res)
  ix=as.numeric(rownames(res))
  pval=res[,'PValue']
  pval=pval[sort(ix,decreasing = F,index.return=T)$ix]
  pval
}


.runDESeq2 <- function(dat) {
  require(DESeq2)
  cond=factor(dat$designs)
  dds=DESeqDataSetFromMatrix(dat$counts, DataFrame(cond), ~ cond)
  if(!dat$sf) sizeFactors(dds)=rep(1,ncol(dat$counts))
  dds=DESeq(dds, quiet=TRUE)
  res=results(dds)
  pval=res$pvalue
  pval[is.na(pval)]=1
  pval
}



.runPois.ztest<-function(dat){
  sfs=colSums(dat$counts);sfs=sfs/min(sfs)
  if(dat$sf) dat$counts=sweep(dat$counts,2,sfs,FUN='/')
  n0=sum(dat$designs==0)
  n1=sum(dat$designs==1)
  m0=rowMeans(dat$counts[,dat$designs==0])
  m1=rowMeans(dat$counts[,dat$designs==1])
  n=nrow(dat$counts)
  pval=rep(1,n)
  for(i in 1:n){
    z=(m1[i]-m0[i])/sqrt(m1[i]/n1+m0[i]/n0)
    pval[i]=2*pnorm(-abs(z))
  }
  pval
}



.runPois.ztest2<-function(dat){
  sfs=colSums(dat$counts);sfs=sfs/min(sfs)
  if(dat$sf) dat$counts=sweep(dat$counts,2,sfs,FUN='/')
  n0=sum(dat$designs==0)
  n1=sum(dat$designs==1)
  m0=rowMeans(dat$counts[,dat$designs==0])
  m1=rowMeans(dat$counts[,dat$designs==1])
  n=nrow(dat$counts)
  pval=rep(1,n)
  for(i in 1:n){
    z=(sqrt(m1[i])-sqrt(m0[i]))/2/sqrt(1/n1+1/n0)
    pval[i]=2*pnorm(-abs(z))
  }
  pval
}




.runPois.glm<-function(dat){
  sfs=colSums(dat$counts)
  n=nrow(dat$counts)
  pval=rep(1,n)
  for(i in 1:n){
    if(dat$sf) fit=glm(dat$counts[i,]~dat$designs,family='poisson',offset=log(sfs))
    else fit=glm(dat$counts[i,]~dat$designs,family='poisson')
    pval[i]=summary(fit)$coefficients[2,4]
  }
  pval
}



.runNB.glm<-function(dat){
  require(MASS)
  n=nrow(dat$counts)
  pval=rep(1,n)
  sfs=colSums(dat$counts)
  for(i in 1:n){
    if(dat$sf) fit=try(glm.nb(dat$counts[i,]~dat$designs,offset=log(sfs)),silent=T)
    else fit=try(glm.nb(dat$counts[i,]~dat$designs),silent=T)
    if(class(fit)!='try-error'){
      pval[i]=summary(fit)$coefficients[2,4]
    }
  }
  pval
}




#######################################
### DE for circRNA (used in realdata)



runedgeR <- function(dat) {
  require(edgeR)
  m=rowMeans(dat$counts)
  id=m>dat$cutoff
  dat$counts=dat$counts[id,]
  dat$gcirc=dat$gcirc[id]
  if(dat$circ.method=='CIRI'){
    dat$nonjuncread=dat$nonjuncread[id,]
    dat$CLR=dat$CLR[id,]
  }
  d=DGEList(counts=dat$counts, group=dat$designs)
  if(dat$sf){
    d=calcNormFactors(d)
  }else{
    d=calcNormFactors(d,method='none')
  }
  d=estimateCommonDisp(d)
  d=estimateTagwiseDisp(d)
  fit=exactTest(d)
  res=topTags(fit, n=nrow(dat$counts))
  res=as.data.frame(res)
  ix=as.numeric(rownames(res))
  pval=res[,'PValue']
  lfc=res[,'logFC']
  fdr=res[,'FDR']
  if(dat$circ.method=='CIRI'){
    res=cbind(as.data.frame(dat$gcirc[ix]),dat$counts[ix,],dat$nonjuncread[ix,],dat$CLR[ix,],lfc=lfc,pval=pval, fdr=fdr)
  }else{
    res=cbind(as.data.frame(dat$gcirc[ix]),dat$counts[ix,],lfc=lfc,pval=pval, fdr=fdr)
  }
  rownames(res)=ix
  res=res[order(res$pval,decreasing = F),]
  res
}




runDESeq2 <- function(dat) {
  require(DESeq2)
  m=rowMeans(dat$counts)
  id=m>dat$cutoff
  dat$counts=dat$counts[id,]
  dat$gcirc=dat$gcirc[id]
  if(dat$circ.method=='CIRI'){
    dat$nonjuncread=dat$nonjuncread[id,]
    dat$CLR=dat$CLR[id,]
  }
  cond=factor(dat$designs)
  dds=DESeqDataSetFromMatrix(dat$counts, DataFrame(cond), ~ cond)
  if(!dat$sf) sizeFactors(dds)=rep(1,ncol(dat$counts))
  dds=DESeq(dds, quiet=TRUE)
  res=results(dds)
  pval=res$pvalue
  pval[is.na(pval)]=1
  fdr=res$padj
  fdr[is.na(fdr)]=1
  lfc=res$log2FoldChange
  if(dat$circ.method=='CIRI'){
    res=cbind(as.data.frame(dat$gcirc),dat$counts,dat$nonjuncread,dat$CLR,lfc=lfc,pval=pval, fdr=fdr)
  }else{
    res=cbind(as.data.frame(dat$gcirc),dat$counts,lfc=lfc,pval=pval, fdr=fdr)
  }
  res=res[order(res$pval,decreasing = F),]
  res
}


runPois.ztest<-function(dat){
  counts=dat$counts
  m=rowMeans(dat$counts)
  id=m>dat$cutoff
  dat$counts=dat$counts[id,]
  dat$gcirc=dat$gcirc[id]
  counts=counts[id,]
  if(dat$circ.method=='CIRI'){
    dat$nonjuncread=dat$nonjuncread[id,]
    dat$CLR=dat$CLR[id,]
  }
  sfs=colSums(dat$counts);sfs=sfs/min(sfs)
  if(dat$sf) dat$counts=sweep(dat$counts,2,sfs,FUN='/')
  n0=sum(dat$designs==0)
  n1=sum(dat$designs==1)
  m0=apply(dat$counts[,dat$designs==0],1,mean)
  m1=apply(dat$counts[,dat$designs==1],1,mean)
  n=nrow(dat$counts)
  pval=rep(1,n)
  for(i in 1:n){
    z=(m1[i]-m0[i])/sqrt(m1[i]/n1+m0[i]/n0)
    pval[i]=2*pnorm(-abs(z))
  }
  fdr=p.adjust(pval,method='fdr')
  lfc=log((m1+1)/(m0+1),2)
  if(dat$circ.method=='CIRI'){
    res=cbind(as.data.frame(dat$gcirc),counts,dat$nonjuncread,dat$CLR,lfc=lfc,pval=pval, fdr=fdr)
  }else{
    res=cbind(as.data.frame(dat$gcirc),counts,lfc=lfc,pval=pval, fdr=fdr)
  }
  res=res[order(res$pval,decreasing = F),]
  res
}



runPois.glm<-function(dat){
  counts=dat$counts
  m=rowMeans(dat$counts)
  id=m>dat$cutoff
  dat$counts=dat$counts[id,]
  dat$gcirc=dat$gcirc[id]
  counts=counts[id,]
  if(dat$circ.method=='CIRI'){
    dat$nonjuncread=dat$nonjuncread[id,]
    dat$CLR=dat$CLR[id,]
  }
  sfs=colSums(dat$counts)
  n=nrow(dat$counts)
  pval=rep(1,n)
  for(i in 1:n){
    if(dat$sf) fit=glm(dat$counts[i,]~dat$designs,family='poisson',offset=log(sfs))
    else fit=glm(dat$counts[i,]~dat$designs,family='poisson')
    pval[i]=summary(fit)$coefficients[2,4]
    #pval.fit[i]=1-pchisq(summary(fit)$deviance,summary(fit)$df.residual)
  }
  
  fdr=p.adjust(pval,method='fdr')
  lfc=log((m1+1)/(m0+1),2)
  if(dat$circ.method=='CIRI'){
    res=cbind(as.data.frame(dat$gcirc),counts,dat$nonjuncread,dat$CLR,lfc=lfc,pval=pval, fdr=fdr)
  }else{
    res=cbind(as.data.frame(dat$gcirc),counts,lfc=lfc,pval=pval, fdr=fdr)
  }
  res=res[order(res$pval,decreasing = F),]
  res
}


runNB.glm<-function(dat){
  counts=dat$counts
  m=rowMeans(dat$counts)
  id=m>dat$cutoff
  dat$counts=dat$counts[id,]
  dat$gcirc=dat$gcirc[id]
  counts=counts[id,]
  if(dat$circ.method=='CIRI'){
    dat$nonjuncread=dat$nonjuncread[id,]
    dat$CLR=dat$CLR[id,]
  }
  sfs=colSums(dat$counts)
  n=nrow(dat$counts)
  pval=rep(1,n)
  for(i in 1:n){
    if(dat$sf) fit=try(glm.nb(dat$counts[i,]~dat$designs,offset=log(sfs)),silent=T)
    else fit=try(glm.nb(dat$counts[i,]~dat$designs),silent=T)
    if(class(fit)!='try-error'){
      pval[i]=summary(fit)$coefficients[2,4]
      #pval.fit[i]=1-pchisq(summary(fit)$deviance,summary(fit)$df.residual)
    }
  }
  fdr=p.adjust(pval,method='fdr')
  lfc=log((m1+1)/(m0+1),2)
  if(dat$circ.method=='CIRI'){
    res=cbind(as.data.frame(dat$gcirc),counts,dat$nonjuncread,dat$CLR,lfc=lfc,pval=pval, fdr=fdr)
  }else{
    res=cbind(as.data.frame(dat$gcirc),counts,lfc=lfc,pval=pval, fdr=fdr)
  }
  res=res[order(res$pval,decreasing = F),]
  res
}










