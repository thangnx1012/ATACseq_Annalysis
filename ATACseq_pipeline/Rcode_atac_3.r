library(ATACseqQC)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPpeakAnno)
library(MotifDb)

samplelabel <- gsub(pattern='.d$',"",list.files(getwd(), pattern='*.d$'))
samplelabel
main_path <- gsub(pattern='/02_clean',"",getwd())
bam_after_shift <- paste(main_path,'/03_shift/',samplelabel,'.aligned.2.clean.shifted.bam',sep='')

bamclean <- paste(main_path,'/02_clean/',samplelabel,'.aligned.2.clean.bam',sep='')
pdf('02_fragment_size_dist_after_QC.pdf')
fragSize <- fragSizeDist(bamclean, bamclean.labels)
dev.off()


motif_interest <- motif_interest <- 'CTCF'

gal3 <- readGAlignments(bam_after_shift)
setwd(paste(main_path,'/99_plot',sep=''))

txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
pt <- PTscore(gal3, txs)

pdf('06_PTscores.pdf')
plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
dev.off()

###Nucleosome Free Regions (NFR) score

nfr <- NFRscore(gal3, txs)
pdf('07_Nuc_free_scores.pdf')
plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()	

###Plot Footprints


CTCF <- query(MotifDb, c(motif_interest))
CTCF <- as.list(CTCF)
seqlev <- c(paste("chr",1:22, sep=''),"chrX","chrY") 

nucleosome_free_bamfile <- paste(main_path,'/04_splited/',samplelabel,'.nucFree.bam',sep='')

pdf(paste(motif_interest,".footprint.pdf",sep=''))
factorFootprints(nucleosome_free_bamfile, pfm=CTCF[[1]], genome=Hsapiens, min.score="90%", seqlev=seqlev, upstream=100, downstream=100)
dev.off()
q()

