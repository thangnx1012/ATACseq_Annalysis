samplelabel <- gsub(pattern='.d$',"",list.files(getwd(), pattern='*.d$'))
samplelabel
main_path <- gsub(pattern='/02_clean',"",getwd())
chrom <- c(paste("chr",1:22, sep=''),"chrX","chrY")

library(ggplot2)
library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)
library(magrittr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(htmlwidgets)
library(rGREAT)
library(Rsamtools)
library(GenomicAlignments)


blkList <- import.bed("/data/rama/labMembers/gg699/atac_pipeline_test/blacklist/ENCFF001TDO.bed.gz")

###Retrieving insert sizes

setwd(paste(main_path,'/03_shift',sep=''))
shiftedBAM <- list.files(path = getwd(),pattern = "\\.shifted.bam$")

atacReads <- readGAlignmentPairs(shiftedBAM, param = ScanBamParam(mapqFilter = 1, 
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE), what = c("qname", 
        "mapq", "isize")))

####Since properly paired reads will have same insert size length we extract insert sizes from read1.

atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)

###Subsetting ATAC-seq reads files by insert sizes

setwd(paste(main_path,'/04_splited',sep=''))

atacReads_Open <- atacReads[insertSizes < 100, ]
atacReads_MonoNuc <- atacReads[insertSizes > 180 & insertSizes < 240, ]
atacReads_diNuc <- atacReads[insertSizes > 315 & insertSizes < 437, ]

export(atacReads_Open, paste(samplelabel,'.nucFree.bam',sep=''), format = "bam")
export(atacReads_MonoNuc,  paste(samplelabel,'.monoNuc.bam',sep=''), format = "bam")
export(atacReads_diNuc, paste(samplelabel,'.diNuc.bam',sep=''), format = "bam")

###Plotting more QC plots

setwd(paste(main_path,'/99_plot',sep=''))

bam_before_QC <- paste(main_path,'/01_bwa/',samplelabel,'.aligned.1.bam',sep='')
bam_after_QC <- paste(main_path,'/02_clean/',samplelabel,'.aligned.2.clean.bam',sep='')
bam_after_shift <- paste(main_path,'/03_shift/',samplelabel,'.aligned.2.clean.shifted.bam',sep='')

pdf('03_Distribution_of_mapped_reads_before_QC.pdf')
ggplot(subset(idxstatsBam(bam_before_QC),idxstatsBam(bam_before_QC)$seqnames==chrom), aes(seqnames, mapped, fill = seqnames)) + 
    geom_bar(stat = "identity") + coord_flip()
dev.off()

pdf('03_Distribution_of_mapped_reads_after_QC.pdf')
ggplot(subset(idxstatsBam(bam_after_QC),idxstatsBam(bam_after_QC)$seqnames==chrom), aes(seqnames, mapped, fill = seqnames)) + 
    geom_bar(stat = "identity") + coord_flip()
dev.off()

pdf('03_Distribution_of_mapped_reads_after_shift.pdf')
ggplot(subset(idxstatsBam(bam_after_shift),idxstatsBam(bam_after_shift)$seqnames==chrom), aes(seqnames, mapped, fill = seqnames)) + 
    geom_bar(stat = "identity") + coord_flip()
dev.off()

fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
    Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
    Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
    geom_line()
	
pdf('04_Insert_sizes_regions_labeled.pdf')
fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315, 
    437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + 
    theme_bw()
	
dev.off()

pdf('04_Insert_sizes_log_regions_labeled.pdf')
fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180, 
    247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") + 
    geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()
dev.off()

##Peak calling with MACS2

setwd(paste(main_path,'/05_peaks',sep=''))

###Un-shifted file

system(paste('macs2 callpeak -t ',bam_after_QC,
' -f BAMPE --nomodel --broad -n ',samplelabel,'.clean_notshifted',sep=''))

###Shifted File
system(paste('macs2 callpeak -t ',bam_after_shift,
' -f BAMPE --nomodel --broad -n ',samplelabel,'.shifted',sep=''))

###Blacklisting peaks 
####Using substract from bedtools all blacklisted regions are removed from peak filel
####and labeled with a suffix  blkList


system(paste('bedtools subtract -a ' , samplelabel,'.clean_notshifted_peaks.bed',
' -b /data/rama/labMembers/gg699/atac_pipeline_test/blacklist/ENCFF001TDO.bed.gz > ',
samplelabel,'.clean_notshifted_peaks_blkList.bed',sep=''))

system(paste('bedtools subtract -a ' , samplelabel,'.clean_notshifted_peaks.bed',
' -b /data/rama/labMembers/gg699/atac_pipeline_test/blacklist/ENCFF001TDO.bed.gz > ',
samplelabel,'.clean_notshifted_peaks_blkList.bed',sep=''))

system(paste('bedtools subtract -a ' , samplelabel,'.shifted_peaks.bed',
' -b ',bam_after_QC,' -counts > ',
samplelabel,'.shifted_peaks_blkList.cts',sep=''))

system(paste('bedtools coverage -a ' , samplelabel,'.clean_notshifted_peaks.bed',
' -b ',bam_after_QC,' -counts > ',
samplelabel,'.clean_notshifted_peaks_blkList.cts',sep=''))

##QC and Annotation of peaks

###For this process there is no need to load the blacklisted peaks
###The script takes care of that, but for other purposes/software you may need them

chrom <- c(paste("chr",1:22, sep=''),"chrX","chrY")

openRegionPeaks_notshifted <- paste(samplelabel,'.clean_notshifted_peaks.bed',sep='')
openRegionPeaks_shifted <- paste(samplelabel,'.shifted_peaks.bed',sep='')

qcRes_notshifted <- ChIPQCsample(bam_after_QC, peaks = openRegionPeaks_notshifted, annotation = "hg19", chromosomes = NULL, blacklist = blkList, verboseT = F)
qcRes_shifted <- ChIPQCsample(bam_after_shift, peaks = openRegionPeaks_shifted, annotation = "hg19", chromosomes = NULL, blacklist = blkList, verboseT = F)


QCmetrics_notshifted <- QCmetrics(qcRes_notshifted) %>% t %>% data.frame %>% dplyr:::select(Reads, starts_with(c("Filt")), 
    starts_with(c("RiP")), starts_with(c("RiBL"))) %>% datatable(rownames = NULL)	

QCmetrics_shifted <- QCmetrics(qcRes_shifted) %>% t %>% data.frame %>% dplyr:::select(Reads, starts_with(c("Filt")), 
    starts_with(c("RiP")), starts_with(c("RiBL"))) %>% datatable(rownames = NULL)	

setwd(paste(main_path,'/99_plot',sep=''))

saveWidget(QCmetrics_notshifted, "QCmetrics_notshifted.html",selfcontained = FALSE)
saveWidget(QCmetrics_shifted, "QCmetrics_shifted.html",selfcontained = FALSE)

flagtagcounts_notshifted <- flagtagcounts(qcRes_notshifted) %>% t %>% data.frame %>% mutate(Dup_Percent = (DuplicateByChIPQC/Mapped) * 
100) %>% dplyr:::select(Mapped, Dup_Percent) %>% datatable(rownames = NULL)

flagtagcounts_shifted <- flagtagcounts(qcRes_shifted) %>% t %>% data.frame %>% mutate(Dup_Percent = (DuplicateByChIPQC/Mapped) * 
100) %>% dplyr:::select(Mapped, Dup_Percent) %>% datatable(rownames = NULL)

saveWidget(flagtagcounts_notshifted,"flagtagcounts_notshifted.html",selfcontained = FALSE)	
saveWidget(flagtagcounts_shifted,"flagtagcounts_shifted.html",selfcontained = FALSE)

qcRes_notshifted_2 <- granges(qcRes_notshifted[seqnames(qcRes_notshifted) %in% chrom])
qcRes_shifted_2 <- granges(qcRes_shifted[seqnames(qcRes_shifted) %in% chrom])

blkList_report <- rbind(data.frame(Blacklisted = sum(qcRes_notshifted_2 %over% blkList), Not_Blacklisted = sum(!qcRes_notshifted_2 %over% 
    blkList)), data.frame(Blacklisted = sum(qcRes_shifted_2 %over% blkList), Not_Blacklisted = sum(!qcRes_shifted_2 %over% 
    blkList)))
rownames(blkList_report) <- c('not_shifted','shifted')

setwd(paste(main_path,'/06_annot',sep=''))

write.table(blkList_report,'Blacklist_report.txt')

qcRes_shifted_2_filtered <- qcRes_shifted_2[!qcRes_shifted_2 %over% blkList]
qcRes_notshifted_2_filtered <- qcRes_notshifted_2[!qcRes_notshifted_2 %over% blkList]


shifted_filteredAnno <- annotatePeak(qcRes_shifted_2_filtered, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
notshifted_filteredAnno <- annotatePeak(qcRes_notshifted_2_filtered, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)

###Print annotation table for shifted and not shifted

write.table(shifted_filteredAnno,"shifted_filteredAnno.txt")
write.table(notshifted_filteredAnno,"notshifted_filteredAnno.txt")


setwd(paste(main_path,'/06_annot',sep=''))

###Annotation plots

pdf('05_plotAnnoPie_shifted.pdf')
plotAnnoPie(shifted_filteredAnno)
dev.off()
pdf('05_plotAnnoBar_shifted.pdf')
plotAnnoBar(shifted_filteredAnno)
dev.off()
pdf('05_upsetplot_shifted.pdf')
upsetplot(shifted_filteredAnno)
dev.off()

pdf('05_plotAnnoPie_notshifted.pdf')
plotAnnoPie(notshifted_filteredAnno)
dev.off()
pdf('05_plotAnnoBar_notshifted.pdf')
plotAnnoBar(notshifted_filteredAnno)
dev.off()
pdf('05_upsetplot_notshifted.pdf')
upsetplot(notshifted_filteredAnno)
dev.off()

###Retrieving annotated Nucleosome-free regions
####Only done with shifted alignment

shiftedGranges_Anno <- as.GRanges(shifted_filteredAnno)
TSS_shiftedGranges_Anno <- shiftedGranges_Anno[abs(shiftedGranges_Anno$distanceToTSS) < 
    500]

####Functional Analysis of Nucleosome-free regions - 1
####This needs to be done in a subset of genes in the genome and not the 
####entire genome because of saturation of hypergeometric tests

seqlevelsStyle(qcRes_shifted_2_filtered) <- "UCSC"

great_Job <- submitGreatJob(qcRes_shifted_2_filtered, species = "hg19")

availableCategories(great_Job)

great_ResultTable_GO = getEnrichmentTables(great_Job, category = "GO")
names(great_ResultTable_GO)
save(great_ResultTable_GO, file = "Great_Results_GO.RData")

great_ResultTable_phenodata_humandisease = getEnrichmentTables(great_Job, category = "Phenotype Data and Human Disease")
names(great_ResultTable_phenodata_humandisease)
save(great_ResultTable_phenodata_humandisease, file = "Great_Results_phenodata_humandisease.RData")

great_ResultTable_pathway = getEnrichmentTables(great_Job, category = "Pathway Data")
names(great_ResultTable_pathway)
save(great_ResultTable_pathway, file = "Great_Results_pathway.RData")

great_ResultTable_expression = getEnrichmentTables(great_Job, category = "Gene Expression")
names(great_ResultTable_expression)
save(great_ResultTable_expression, file = "Great_Results_expression.RData")

great_ResultTable_regmotifs = getEnrichmentTables(great_Job, category = "Regulatory Motifs")
names(great_ResultTable_regmotifs)
save(great_ResultTable_regmotifs, file = "Great_Results_regmotifs.RData")

great_ResultTable_genefams = getEnrichmentTables(great_Job, category = "Gene Families")
names(great_ResultTable_genefams)
save(great_ResultTable_genefams, file = "Great_Results_genefams.RData")

###Make bedfiles just with peaks at promoter regions

####This peak files are useful to run with HOMER


shift_Anno2 <- read.table("shifted_filteredAnno.txt", head=T)
promoter_regions <- c('Promoter (<=1kb)','Promoter (1-2kb)','Promoter (2-3kb)')
t1 <- subset(shift_Anno2,  shift_Anno2$annotation == promoter_regions[1])
t2 <- subset(shift_Anno2,  shift_Anno2$annotation == promoter_regions[2])
t3 <- subset(shift_Anno2,  shift_Anno2$annotation == promoter_regions[3])
shift_Anno2_promoter<- rbind(t1,t2,t3)
rm(t1)
rm(t2)
rm(t3)
write.table(shift_Anno2_promoter,"filteredAnno_promoter.txt")

setwd(paste(main_path,'/05_peaks',sep=''))

bed_promoters <- shift_Anno2_promoter[,c('seqnames','start','end')]
colnames(bed_promoters) <- c('chrom','chromStart','chromEnd')
write.table(bed_promoters,"bed_promoters.txt", row.names=F, col.names=F)
q()