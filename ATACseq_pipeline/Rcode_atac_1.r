library(ATACseqQC)
samplelabel <- gsub(pattern='.d',"",list.files(getwd(), pattern='*.d$'))
samplelabel
main_path <- gsub(pattern='/02_clean',"",getwd())

setwd(paste(main_path,'/01_bwa/',sep=''))

sortedBAM <- list.files(path = paste(main_path,'/01_bwa/',sep=''),pattern = "\\.aligned.1.bam$")
sortedBAM.labels <- sub(".bam", "", basename(sortedBAM))
rawBAM <- list.files(path = paste(main_path,'/01_bwa/',sep=''),pattern = "\\.aligned.0.bam$")
rawBAM.labels <- sub(".bam", "", basename(rawBAM))

##Quality control on Bam file: Removing duplicates...etc

bamQC(sortedBAM, index = sortedBAM, mitochondria = "chrM",
outPath = sub(".bam", ".clean.bam", basename(sortedBAM)),
doubleCheckDup = TRUE)

##Plots for Quality control
###Estimate library complexity and compare noM with clean files

pdf('01_lib_complexity_raw.pdf')
estimateLibComplexity(readsDupFreq(rawBAM))
dev.off()

pdf('01_lib_complexity_before_QC.pdf')
estimateLibComplexity(readsDupFreq(sortedBAM))
dev.off()

pdf('01_lib_complexity_after_QC.pdf')
estimateLibComplexity(readsDupFreq(sub(".bam", ".clean.bam", basename(sortedBAM))))
dev.off()

bamfile <- sub(".bam", ".clean.bam", basename(sortedBAM))
bamfile.labels <- sub(".bam", "", basename(bamfile))

###Fragment size distribution

pdf('02_fragment_size_dist_raw.pdf')
fragSize <- fragSizeDist(rawBAM, rawBAM.labels)
dev.off()

pdf('02_fragment_size_dist_before_QC.pdf')
fragSize <- fragSizeDist(sortedBAM, sortedBAM.labels)
dev.off()

pdf('02_fragment_size_dist_after_QC.pdf')
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()

###Move plots to plot folder

system(paste('mv *.pdf ',main_path,'/99_plot/.',sep=''))

###Move clean bam to folder

system(paste('mv ',bamfile,' ',main_path,'/02_clean/.',sep=''))
system(paste('mv ',bamfile,'.bai ',main_path,'/02_clean/.',sep=''))
q()