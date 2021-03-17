source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
ls("package:BSgenome.Hsapiens.UCSC.hg19")

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
