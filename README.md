# circMeta
A unified computational framework for genomic feature annotation, differential expression analysis of circular RNAs

## Description
circMeta is a unified computational framework for circRNA analyses. circMeta mainly includes three function modules: (i) provide a comprehensive genomic feature annotation related to circRNA biogenesis, including length of introns flanking circularized exons, repetitive elements such as Alu and SINEs, competition score for forming circulation and RNA editing in back-spliced flanking introns (ii) develop a two-stage DE approach of circRNAs based on splicing junction reads (iii) develop a Bayesian hierarchical model for DE analysis of circRNAs based the ratio of circular reads to linear reads in spliced sites.
circMeta mainly consists of four modules: circClass, circFeature, circJuncDE, circCLRDE.

# Maintainer
Li Chen <li.chen@auburn.edu>


# Install circMeta
```r
install.packages("devtools")
library(devtools)
install_github("lichen-lab/circMeta")
```


# Descriptions for circJuncDE

## Usage
circJuncDE(files,designs, circ.method=c('findcirc','CIRCexplorer','CIRI'),
                 DE.method=c('pois.ztest','DESeq2','edgeR','pois.glm','nb.glm'),
                 gene=NULL,gexon=NULL,
                 cutoff=2,sf=TRUE)

## Arguments
*  files: circRNA output files from 'findcirc', 'CIRCexplorer' or 'CIRI'
*  designs: design matrix for files. For example, two group design is c(0,0,1,1)
*  circ.method: one of 'findcirc','CIRCexplorer','CIRI' should be specified
*  DE.method: one of 'DESeq2','edgeR','pois.glm','pois.ztest','nb.glm'. Default is 'pois.ztest'
*  gene: gene annotation is optional. Default is NULL
*  gexon: exon annotation is optional. Default is NULL
*  cutoff: cutoff for circRNAs. Default is 2
*  sf: sequencing depth is adjusted for all juntion reads. Default is TRUE
                 
                
## Output values
* a GRange object contains genomic coordinates of circRNAs along with junction reads, fold change, p-values, FDRs.


# Descriptions for circCLEDE

## Usage
circCLEDE(x1,n1,x2,n2,DE.method=c('wald','lr','fisher','chisq'),is.shrink=TRUE,is.equalrho=FALSE,is.peudo=TRUE)

## Arguments
*  x1: junction reads for circRNAs in condition1
*  n1: junction and linear reads for spliced sites in condition1
*  x2: junction reads for circRNAs in condition2
*  n2: junction and linear reads for spliced sites in condition2
*  DE.method: one of 'wald','lr','fisher','chisq'
*  is.shrink: use shrunk rho (dispersion) estimate for 'wald' or 'lr'. Default is 'wald'
*  is.equalrho: rho is the same for two conditions. Default is FALSE
*  is.peudo: add peudo count 1 to circRNAs where either n1 or n2 is 0

## Output values
* a GRange object contains genomic coordinates of circRNAs along with junction reads, nonjunction reads, CLR, p-values, FDRs.


# Descriptions for circClass

## Usage
circClass(files,circ.method=c('findcirc','CIRCexplorer','CIRI'),gene=NULL,gexon=NULL)

## Arguments
*  files: circRNA output files from 'findcirc', 'CIRCexplorer' or 'CIRI'
*  circ.method: one of 'findcirc','CIRCexplorer','CIRI' should be specified
*  gene: gene annotation is optional. Default is NULL
*  gexon: exon annotation is optional. Default is NULL

## Output values
* a GRange object contains genomic coordinates of circRNAs along with host gene annotation.



# Descriptions for circFeature

## Usage
circFeature(files,circ.method=c('findcirc','CIRCexplorer','CIRI'),
                      gene,gexon,gintron,galu,gsine,
                      atoi=NULL,DARNED=NULL,RADAR=NULL)

## Arguments
*  files: circRNA output files from 'findcirc', 'CIRCexplorer' or 'CIRI'
*  circ.method: one of 'findcirc','CIRCexplorer','CIRI' should be specified
*  gene: genes of hg19 or mm9 in GRange object 
*  gexon: exon of hg19 or mm9 in GRange object 
*  gintron: intron of hg19 or mm9 in GRange object 
*  galu: Alu of hg19 or mm9 in GRange object 
*  gsine: SINE of hg19 or mm9 in GRange object 
*  atoi: A-to-I editing in GRange object. Default is NULL
*  DARNED: DARNED RNA editing in GRange object. Default is NULL
*  RADAR: RADAR A-to-I editing in GRange object. Default is NULL 


## Output values
* a GRange object contains genomic coordinates of circRNAs along with all genomic features listed as the parameters



# Descriptions for circLinear

## Usage
circLinear(gcirc,bamfiles)


## Arguments
*  gcirc: a GRange object contains genomic coordinates of circRNAs 
*  bamfiles: bam files containing mapped linear reads e.g. tophat

## Output values
* a GRange object contains genomic coordinates of circRNAs along with junction reads and nonjunction reads in back-spliced sites





# Examples

Files used in the examples could be downloaded from [link](https://drive.google.com/open?id=1Eaxbx8w33_Mjxcd9KPIsrSkflZSuAIHV)

Use circJuncDE to detect differential expressed circRNAs based on back-spliced junction reads
```r
library(circMeta)
CIRI.files=c(
  "data/frontcortex/CIRI/CIRI1.circ",
  "data/frontcortex/CIRI/CIRI2.circ",
  "data/cerebellum/CIRI/CIRI1.circ",
  "data/cerebellum/CIRI/CIRI2.circ",
  "data/diencephalon/CIRI/CIRI1.circ",
  "data/diencephalon/CIRI/CIRI2.circ"
)
designs=c(0,0,1,1)
fdr.level=0.05
res1=circJuncDE(CIRI.files[c(1,2,3,4)],designs,circ.method='CIRI',DE.method='pois.ztest')
id1=rownames(res1)[res1$fdr<fdr.level]

res2=circJuncDE(CIRI.files[c(1,2,3,4)],designs,circ.method='CIRI',DE.method='edgeR')
id2=rownames(res2)[res2$fdr<fdr.level]

res3=circJuncDE(CIRI.files[c(1,2,3,4)],designs,circ.method='CIRI',DE.method='DESeq2')
id3=rownames(res3)[res3$fdr<fdr.level]


id1=rownames(res1)[res1$fdr<fdr.level]
length(id1)
id2=rownames(res2)[res2$fdr<fdr.level]
length(id2)
id3=rownames(res3)[res3$fdr<fdr.level]
length(id3)
length(intersect(id1,id3))
length(intersect(id2,id3))
length(intersect(id1,id2))

```

Use circCLRDE to detect differential expressed circRNAs based on Circular-to-Linear Ratio (CLR)
```r
res=circJuncDE(CIRI.files[c(1,2,3,4)],designs,circ.method='CIRI',DE.method='pois.ztest')
tmp=getCLR(res)
x1=tmp$x1
n1=tmp$n1
x2=tmp$x2
n2=tmp$n2

res1=circCLRDE(x1,n1,x2,n2,DE.method='wald',is.shrink=F)
sum(res1$fdr<fdr.level)
res2=circCLRDE(x1,n1,x2,n2,DE.method='wald',is.shrink=T)
sum(res2$fdr<fdr.level)
res3=circCLRDE(x1,n1,x2,n2,DE.method='lr',is.shrink=F)
sum(res3$fdr<fdr.level)
res4=circCLRDE(x1,n1,x2,n2,DE.method='lr',is.shrink=T)
sum(res4$fdr<fdr.level)
```

Use circClass to annotate host genes for circRNAs
```r
load('hg19/Exon.hg19.rda')
load('hg19/Gene.hg19.rda')
load('hg19/Intron.hg19.rda')
gcirc=circClass(files=CIRI.files[1:2],circ.method='CIRI',gene=gene, gexon=gexon)
gcirc
```


Use circFeature to annotate genomic features for circRNAs
```r
load('hg19/Alu.hg19.rda')
load('hg19/SINE.hg19.rda')

# Add Alu/SINE enrichment, Alu/SINE competing score
gcirc=circFeature(files=CIRI.files[1:2],circ.method='CIRI',
                     gene,gexon,gintron,galu,gsine)

# Add A-to-I editing
load('hg19/DARNED.rda')
load('hg19/RADAR.rda')
atoi=getAtoI('hg19/frontalcortex.vcf2')
gcirc=circFeature(files=CIRI.files[1:2],circ.method='CIRI',
                     gene,gexon,gintron,galu,gsine,atoi,DARNED,RADAR)
```

Use circLinear to add linear reads to back-spliced junctions
```r
llinear.files=c(
  'hg19/accepted_hits1.bam',
  'hg19/accepted_hits2.bam'
)
gcirc=circLinear(gcirc,linear.files)
```











