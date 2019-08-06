##################
### A-to-I editing

module load samtools/1.9
module load bcftools/1.9

genome=/home/Shared/JinLab/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

bam1=/home/lchen/ncbi/public/sra/ENCODE/frontalcortex/ENCFF001RN/circexplorer/tophat/accepted_hits.bam
bam2=/home/lchen/ncbi/public/sra/ENCODE/frontalcortex/ENCFF001RO/circexplorer/tophat/accepted_hits.bam

bcftools mpileup -f $genome $bam1 $bam2 | bcftools call -mv -Ob -o frontalcortex.bcf
bcftools view frontalcortex.bcf | bcftools view > frontalcortex.vcf
cut -f1,2,4,5 frontalcortex.vcf > frontalcortex.vcf2

############
### findcirc

fastq=/home/lchen/ncbi/public/sra/ENCODE/diencephalon/ENCFF001RN/ENCFF001RM.fastq.gz
bowtie2Index=/home/Shared/JinLab/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
genome=/home/Shared/JinLab/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
findcirc=/home/lchen/software/find_circ/
bowtie2 -p16 --very-sensitive --score-min=C,-15,0 --mm \
-x $bowtie2Index  -q -U $fastq 2> bowtie2.log  \
| samtools view -hbuS - | samtools sort - bowtie2

samtools view -hf 4 bowtie2.bam | samtools view -Sb - > unmapped.bam


$findcirc"unmapped2anchors.py" unmapped.bam >anchors.fastq
gzip anchors.fastq


mkdir -p circ
bowtie2 -p 16 --score-min=C,-15,0 --reorder --mm \
-q -U anchors.fastq.gz -x $bowtie2Index |\
$findcirc"find_circ.py" \
--genome=$genome \
--prefix= \
--name=\
--stats=circ/stats.txt \
--reads=circ/spliced_reads.fa \
> circ/splice_sites.bed



grep CIRCULAR circ/splice_sites.bed | \
awk '$5>=2' | \
grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | \
$findcirc"maxlength.py" 100000 \
> circ/circ_candidates.bed



################
## circexplorer

fastq=/home/lchen/ncbi/public/sra/ENCODE/diencephalon/ENCFF001RN/ENCFF001RM.fastq.gz
genegtf=/home/Shared/PengLab/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf 
bowtie2Index=/home/Shared/PengLab/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
bowtieIndex=/home/Shared/PengLab/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome
genome=/home/Shared/JinLab/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
CIRCexplorer=/home/lchen/software/CIRCexplorer/circ/CIRCexplorer.py
ref=/home/lchen/projects/hg19/hg19.ref.txt

mkdir tophat
tophat2 -p 15 -n 2 -G $genegtf -o tophat $bowtie2Index  $fastq

bamToFastq -i tophat/unmapped.bam -fq tophat/unmapped.fastq

tophat2 -o tophat_fusion -p 15 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search $bowtieIndex tophat/unmapped.fastq
$CIRCexplorer -f tophat_fusion/accepted_hits.bam -g $genome -r $ref

################
## CIRI

fastq=/home/lchen/ncbi/public/sra/ENCODE/diencephalon/ENCFF001RN/ENCFF001RM.fastq.gz
name=ENCFF001RO
in=$name.sam
out=$name.circ
genome=/home/Shared/JinLab/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
gtf=/home/Shared/JinLab/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
circ=/home/lchen/software/CIRI/CIRI2.pl

bwa mem $genome $fastq -T 19 > $in
perl $circ -I $in -O $out -F  $genome -A $gtf





































