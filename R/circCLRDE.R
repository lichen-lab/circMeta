circCLRDE<-function(x1,n1,x2,n2,
                    DE.method=c('wald','lr','fisher','chisq'),
                    is.shrink=T,is.equalrho=F,is.peudo=T) {
  
  DE.method=match.arg(DE.method)
  nrep1=ncol(x1)
  nrep2=ncol(x2)
  
  if(is.peudo){
    ix=rowSums(n1==0)>0
    n1[ix,]=n1[ix,]+1
    ix=rowSums(n2==0)>0
    n2[ix,]=n2[ix,]+1
  }
  
  if(DE.method=='chisq'){
    res=chisqTest(x1,n1,x2,n2)
    return(res)
  }else if(DE.method=='fisher'){
    res=fisherTest(x1,n1,x2,n2)
    return(res)
  }
  
  message('estimate means......')
  est.mu1=estMu(x1, n1)
  est.mu2=estMu(x2, n2)
  est.mu=estMu(x1+x2, n1+n2)
  
  message('estimate rhos......')
  if(is.equalrho | nrep1==1 | nrep2==1) {
    rho1=rho2=estRhoBayes(cbind(x1,x2), cbind(n1,n2), cbind(est.mu1, est.mu2))
  }else {
    if(is.shrink){
      rho1=estRhoBayes(x1, n1, est.mu1)
      rho2=estRhoBayes(x2, n2, est.mu2)
    }else{
      rho1=estRho(x1, n1)
      rho2=estRho(x2, n2)
    }
  }
  
  message('perform statistical test......')
  if(DE.method=='wald'){
    wald=waldTest(n1, est.mu1[,1], rho1, n2, est.mu2[,1], rho2)
    stat=wald$stat; pval=wald$pval; fdr=wald$fdr
  }else if(DE.method=='lr'){
    lr=lrTest(x1, n1, est.mu1, rho1, x2, n2, est.mu2, rho2,est.mu)
    stat=lr$stat; pval=lr$pval; fdr=lr$fdr
  }
  dif=rowMeans(est.mu1-est.mu2)
  
  data.frame(mu1=est.mu1[,1], mu2=est.mu2[,1],dif=dif,rho1=rho1, rho2=rho2,
             stat=stat,pval=pval,fdr=fdr)
}



