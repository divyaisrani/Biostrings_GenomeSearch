# Biostrings_GenomeSearch
search for sequences in genome

Adapted from:
https://bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/GenomeSearching.pdf
https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html

Load the Biostrings human genome package into R
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
ls("package:BSgenome.Hsapiens.UCSC.hg19")

Set genome and all chromosomes as objects
genome <- BSgenome.Hsapiens.UCSC.hg19
genome
chr1 <- genome$chr1
chr2 <- genome$chr2
chr3 <- genome$chr3
chr4 <- genome$chr4
chr5 <- genome$chr5
chr6 <- genome$chr6
chr7 <- genome$chr7
chr8 <- genome$chr8
chr9 <- genome$chr9
chr10 <- genome$chr10
chr11 <- genome$chr11
chr12 <- genome$chr12
chr13 <- genome$chr13
chr14 <- genome$chr14
chr15 <- genome$chr15
chr16 <- genome$chr16
chr17 <- genome$chr17
chr18 <- genome$chr18
chr19 <- genome$chr19
chr20 <- genome$chr20
chr21 <- genome$chr21
chr22 <- genome$chr22
chrX <- genome$chrX
chrY <- genome$chrY
chrM <- genome$chrM

----------------------------------------------
hgenomesearch.R
Use countPattern to return total number of genomic search counts for each sequence in a matrix. You can
run multiple sequences in the same matrix simultaneously, separated by a comma. Remember to also run
each sequence’s reverse complement and add this total to the forward sequence total (genome &amp;
chromosome objects only contain forward strand sequence).

seqmatrix18fwd <- c('GTCGAAGTCGAAGTCGAC','GACGACGTTACGGACGTA','GACGTCGAAGTAGCCGTA','GACGACGCCGATGTAGAA','GAAGCAGTCGACGCCGAA','GACGACGCGGTCTAAGAA','GACGAGGTCGCATAAGTA','GACGCAGTATAGGTCGAA','GACGCAGTATAGGACGAC','GGCGTAGCCGATGTCGCG','GTCGTTGCGGTAGTCGAA')
seqmatrix18rev <- c('GTCGACTTCGACTTCGAC','TACGTCCGTAACGTCGTC','TACGGCTACTTCGACGTC','TTCTACATCGGCGTCGTC','TTCGGCGTCGACTGCTTC','TTCTTAGACCGCGTCGTC','TACTTATGCGACCTCGTC','TTCGACCTATACTGCGTC','GTCGTCCTATACTGCGTC','CGCGACATCGGCTACGCC','TTCGACTACCGCAACGAC')

hgenomesearch <- function(seqmatrix){
for (seq in seqmatrix){
count=countPattern(seq,chr1)+countPattern(seq,chr2)+countPattern(seq,chr3)+countPattern(seq,chr4)+countPattern(seq,chr5)+countPattern(seq,chr6)+countPattern(seq,chr7)+countPattern(seq,chr8)+countPattern(seq,chr9)+countPattern(seq,chr10)+countPattern(seq,chr11)+countPattern(seq,chr12)+countPattern(seq,chr13)+countPattern(seq,chr14)+countPattern(seq,chr15)+countPattern(seq,chr16)+countPattern(seq,chr17)+countPattern(seq,chr18)+countPattern(seq,chr19)+countPattern(seq,chr20)+countPattern(seq,chr21)+countPattern(seq,chr22)+countPattern(seq,chrX)+countPattern(seq,chrY)+countPattern(seq,chrM)
print(count)
}
}

hgenomesearch(seqmatrix18fwd)
hgenomesearch(seqmatrix18rev)

----------------------------------------------
hgenomesearchmismatch.R
Use countPattern with max.mismatch to return total number of genomic search counts with mismatch
tolerance for each sequence in a matrix. You can run multiple sequences in the same matrix simultaneously,
separated by a comma. Remember to also run each sequence’s reverse complement and add this total to the
forward sequence total (genome &amp; chromosome objects only contain forward strand sequence).

seqmatrix18fwd <- c('GTCGAAGTCGAAGTCGAC','GACGACGTTACGGACGTA','GACGTCGAAGTAGCCGTA','GACGACGCCGATGTAGAA','GAAGCAGTCGACGCCGAA','GACGACGCGGTCTAAGAA','GACGAGGTCGCATAAGTA','GACGCAGTATAGGTCGAA','GACGCAGTATAGGACGAC','GGCGTAGCCGATGTCGCG','GTCGTTGCGGTAGTCGAA')
seqmatrix18rev <- c('GTCGACTTCGACTTCGAC','TACGTCCGTAACGTCGTC','TACGGCTACTTCGACGTC','TTCTACATCGGCGTCGTC','TTCGGCGTCGACTGCTTC','TTCTTAGACCGCGTCGTC','TACTTATGCGACCTCGTC','TTCGACCTATACTGCGTC','GTCGTCCTATACTGCGTC','CGCGACATCGGCTACGCC','TTCGACTACCGCAACGAC')

hgenomesearchmismatch <- function(seqmatrix){
for (seq in seqmatrix){
count=countPattern(seq,chr1,max.mismatch=4)+countPattern(seq,chr2,max.mismatch=4)+countPattern(seq,chr3,max.mismatch=4)+countPattern(seq,chr4,max.mismatch=4)+countPattern(seq,chr5,max.mismatch=4)+countPattern(seq,chr6,max.mismatch=4)+countPattern(seq,chr7,max.mismatch=4)+countPattern(seq,chr8,max.mismatch=4)+countPattern(seq,chr9,max.mismatch=4)+countPattern(seq,chr10,max.mismatch=4)+countPattern(seq,chr11,max.mismatch=4)+countPattern(seq,chr12,max.mismatch=4)+countPattern(seq,chr13,max.mismatch=4)+countPattern(seq,chr14,max.mismatch=4)+countPattern(seq,chr15,max.mismatch=4)+countPattern(seq,chr16,max.mismatch=4)+countPattern(seq,chr17,max.mismatch=4)+countPattern(seq,chr18,max.mismatch=4)+countPattern(seq,chr19,max.mismatch=4)+countPattern(seq,chr20,max.mismatch=4)+countPattern(seq,chr21,max.mismatch=4)+countPattern(seq,chr22,max.mismatch=4)+countPattern(seq,chrX,max.mismatch=4)+countPattern(seq,chrY,max.mismatch=4)+countPattern(seq,chrM,max.mismatch=4)
print(count)
}
}

hgenomesearchmismatch(seqmatrix18fwd)
hgenomesearchmismatch(seqmatrix18rev)

----------------------------------------------
hgenomematch.R
Use matchPattern to return information on genomic location for each sequence in a matrix. Use
matchPattern with max.mismatch to find the sequence and genomic locations of mismatches as well. You
can run multiple sequences in the same matrix simultaneously, separated by a comma. Remember to also
run each sequence’s reverse complement (genome &amp; chromosome objects only contain forward strand
sequence).

seqmatrix18fwd <- c('GTCGAAGTCGAAGTCGAC','GACGACGTTACGGACGTA','GACGTCGAAGTAGCCGTA','GACGACGCCGATGTAGAA','GAAGCAGTCGACGCCGAA','GACGACGCGGTCTAAGAA','GACGAGGTCGCATAAGTA','GACGCAGTATAGGTCGAA','GACGCAGTATAGGACGAC','GGCGTAGCCGATGTCGCG','GTCGTTGCGGTAGTCGAA')
seqmatrix18rev <- c('GTCGACTTCGACTTCGAC','TACGTCCGTAACGTCGTC','TACGGCTACTTCGACGTC','TTCTACATCGGCGTCGTC','TTCGGCGTCGACTGCTTC','TTCTTAGACCGCGTCGTC','TACTTATGCGACCTCGTC','TTCGACCTATACTGCGTC','GTCGTCCTATACTGCGTC','CGCGACATCGGCTACGCC','TTCGACTACCGCAACGAC')

hgenomematch <- function(seqmatrix){
for (seq in seqmatrix){
match1=matchPattern(seq,chr1,max.mismatch=2)
match2=matchPattern(seq,chr2,max.mismatch=2)
match3=matchPattern(seq,chr3,max.mismatch=2)
match4=matchPattern(seq,chr4,max.mismatch=2)
match5=matchPattern(seq,chr5,max.mismatch=2)
match6=matchPattern(seq,chr6,max.mismatch=2)
match7=matchPattern(seq,chr7,max.mismatch=2)
match8=matchPattern(seq,chr8,max.mismatch=2)
match9=matchPattern(seq,chr9,max.mismatch=2)
match10=matchPattern(seq,chr10,max.mismatch=2)
match11=matchPattern(seq,chr11,max.mismatch=2)
match12=matchPattern(seq,chr12,max.mismatch=2)
match13=matchPattern(seq,chr13,max.mismatch=2)
match14=matchPattern(seq,chr14,max.mismatch=2)
match15=matchPattern(seq,chr15,max.mismatch=2)
match16=matchPattern(seq,chr16,max.mismatch=2)
match17=matchPattern(seq,chr17,max.mismatch=2)
match18=matchPattern(seq,chr18,max.mismatch=2)
match19=matchPattern(seq,chr19,max.mismatch=2)
match20=matchPattern(seq,chr20,max.mismatch=2)
match21=matchPattern(seq,chr21,max.mismatch=2)
match22=matchPattern(seq,chr22,max.mismatch=2)
matchX=matchPattern(seq,chrX,max.mismatch=2)
matchY=matchPattern(seq,chrY,max.mismatch=2)
matchM=matchPattern(seq,chrM,max.mismatch=2)
print(match1)
print(match2)
print(match3)
print(match4)
print(match5)
print(match6)
print(match7)
print(match8)
print(match9)
print(match10)
print(match11)
print(match12)
print(match13)
print(match14)
print(match15)
print(match16)
print(match17)
print(match18)
print(match19)
print(match20)
print(match21)
print(match22)
print(matchX)
print(matchY)
print(matchM)
}
}

hgenomematch(seqmatrix18fwd)
hgenomematch(seqmatrix18rev)
