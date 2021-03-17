source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
ls("package:BSgenome.Hsapiens.UCSC.hg19&quot")

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

hgenomesearch <- function(seqmatrix){
for (seq in seqmatrix){
count=countPattern(seq,chr1)+countPattern(seq,chr2)+countPattern(seq,chr3)+countPattern(seq,chr4)+countPattern(seq,chr5)+countPattern(seq,chr6)+countPattern(seq,chr7)+countPattern(seq,chr8)+countPattern(seq,chr9)+countPattern(seq,chr10)+countPattern(seq,chr11)+countPattern(seq,chr12)+countPattern(seq,chr13)+countPattern(seq,chr14)+countPattern(seq,chr15)+countPattern(seq,chr16)+countPattern(seq,chr17)+countPattern(seq,chr18)+countPattern(seq,chr19)+countPattern(seq,chr20)+countPattern(seq,chr21)+countPattern(seq,chr22)+countPattern(seq,chrX)+countPattern(seq,chrY)+countPattern(seq,chrM)
print(count)
}
}

hgenomesearch(seqmatrix18fwd)
hgenomesearch(seqmatrix18rev)
