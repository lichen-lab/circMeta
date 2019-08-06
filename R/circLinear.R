circLinear<-function(gcirc,bamfiles){
  gcirc1=gcirc2=gcirc
  end(gcirc1)=start(gcirc1)
  start(gcirc2)=end(gcirc2)
  res1=res2 = matrix(0, nrow = length(gcirc1), ncol = length(bamfiles))
  for (i in 1:length(bamfiles)) {
    reads.ranges = read.BAM(bamfiles[i])
    res1[, i] = countOverlaps(gcirc1, reads.ranges)
    res2[, i] = countOverlaps(gcirc2, reads.ranges)
    rm(reads.ranges)
  }
  res=res1+res2
  colnames(res)= paste('nonjuncread:',1:length(bamfiles),sep='')
  mcols(gcirc)=data.frame(mcols(gcirc),nonjuncread=res)
  gcirc
}

