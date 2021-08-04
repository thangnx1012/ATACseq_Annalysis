
###Specify outbase as an identifier for your sample 
outBase=sample_id


###Specify input and output directories as variables on bash
###Input directory must contain pair-end reads in fq.gz format 

outDir=/full/output/path/$outBase #Modify according to environment
readDir=/full/input/path #Modify according to environment

read1=$readDir/*1.fq.gz
read2=$readDir/*2.fq.gz

conda activate atac_chip_pipeline; ##yml file to replicate environment available in repository


##Defining folder organization

###Reference genome path 
ref=$dbDir/${genome}/${genome}_ERCC92_RGC/${genome}_ERCC92_RGC.fa #Modify accordingly

###Folders
###The pipeline creates a folder structure for each of the results and intermidiate files created 

fastqDir=$outDir/00_fastq
bwaDir=$outDir/01_bwa
bamcleanDir=$outDir/02_clean
bamshiftDir=$outDir/03_shift
bamsplitDir=$outDir/04_splited
peakDir=$outDir/05_peaks
annotation=$outDir/06_annot
bwigDir=$outDir/07_bigwig
plotsDir=$outDir/99_plot
tmpDir=$outDir/02_clean/02_tmp/$outBase
binsDir=$outDir/08_bin

mkdir -p $outDir
mkdir -p $fastqDir
mkdir -p $bwaDir
mkdir -p $bamcleanDir
mkdir -p $peakDir
mkdir -p $annotation
mkdir -p $plotsDir
mkdir -p $bwaDir/tmp
mkdir -p $bwigDir
mkdir -p $tmpDir
mkdir -p $binsDir
mkdir -p $bamshiftDir
mkdir -p $bamsplitDir

##Triming

###FASTQC GZ Trimmed files

##Triming

###Unziping
cd $fastqDir
gunzip $read1 -c > $fastqDir/read1.fq
gunzip $read2 -c > $fastqDir/read2.fq

fqs=$(find $fastqDir -type f -name "*.fq")

###Triming with trim galore
trim_galore --fastqc --paired $fqs

###Removing the un-trimmed reads
rm read1.fq
rm read2.fq

##Aligning
fqcs=$(find $fastqDir -type f -name "*.fq")
bwa mem -R '@RG\tID:1\tSM:'${outBase}'\tPL:illumina\tLB:lib1\tPU:unit1' $ref $fqcs | samtools view -bS -F 0x900 - > $bwaDir/$outBase.aligned.0.bam

samtools view -c -F 0x4 $bwaDir/$outBase.aligned.0.bam > \
    $bwaDir/$outBase.numAligned.txt

samtools sort -T $bwaDir/tmp/ -o $bwaDir/$outBase.aligned.0.bam $bwaDir/$outBase.aligned.0.bam
samtools index $bwaDir/$outBase.aligned.0.bam
touch $bwaDir/$outBase.done

##Cleaning

##Filter BAM file and remove un-scaffolded regions and chrM

samtools view -b -L /data/rama/labMembers/gg699/human_hg19.bed $bwaDir/$outBase.aligned.0.bam > $bwaDir/$outBase.aligned.1.bam
samtools sort -T $bwaDir/tmp/aln.sorted -o $bwaDir/$outBase.aligned.1.bam $bwaDir/$outBase.aligned.1.bam
samtools index $bwaDir/$outBase.aligned.1.bam


cd $bamcleanDir

touch $outBase.d
conda deactivate;
conda activate r_3_6_env; ##yml file to replicate environment available in repository
R CMD BATCH /data/rama/labMembers/gg699/atac_pipeline_test/march/Rcode_atac_1.r;
conda deactivate;

## Filter low-quality reads 
## This step is implemented using a custom tool called FilterBamv2 developed by the Lawrence Lab hosted
## in the Broad Institute server 

conda activate atac_chip_pipeline;

clspth=/data/rama/lawrence/jars/htsjdk-2.0.1.jar
clspth=$clspth:/data/rama/lawrence/cga/trunk/tools/classes

java -Xmx15g -Djava.io.tmpdir=$tmpDir \
    -classpath $clspth org.broadinstitute.cga.tools.seq.FilterBAMv2 \
    -nm 2 -asr 0.67 -mq 30 -dup 0 \
    $bamcleanDir/$outBase.aligned.1.clean.bam \
    $bamcleanDir/$outBase.aligned.2.clean.bam

samtools index $outBase.aligned.2.clean.bam;

#Binning

cd $binsDir

conda activate atac_chip_pipeline;

hg19_2kb_bins=/data/rama/labMembers/gg699/hg19.2kb.windows.bed

bedtools coverage -a $hg19_2kb_bins -b $bamcleanDir/$outBase.aligned.2.clean.bam  -mean > $binsDir/$outBase.2kb.bedgraph
bedtools coverage -a $hg19_2kb_bins -b $bamcleanDir/$outBase.aligned.2.clean.bam  > $binsDir/$outBase.2kb.counts

hg19_200b_bins=/data/rama/labMembers/gg699/hg19.200bp.bed

bedtools coverage -a $hg19_200b_bins -b $bamcleanDir/$outBase.aligned.2.clean.bam  -mean > $binsDir/$outBase.200b.bedgraph
bedtools coverage -a $hg19_200b_bins -b $bamcleanDir/$outBase.aligned.2.clean.bam  > $binsDir/$outBase.200b.counts

##BigWigs created with deeptools

bamCoverage -b $bamcleanDir/$outBase.aligned.2.clean.bam --normalizeUsing CPM --outFileFormat bigwig --binSize 1 -o $bwigDir/$outBase.1.CPM.bw

##Shifting BAM

###This step is crucial for footprinting analysis according to the literature, but it can be skiped if you are not doing that
###The shifting is done by a program within deeptools called alignment sieve that is installed in my folder

alignmentSieve -b $bamcleanDir/$outBase.aligned.2.clean.bam \
-o $bamshiftDir/$outBase.aligned.2.clean.shifted.bam --ATACshift;


###Sorting and indexing new shifted file
cd $bamshiftDir
mkdir tmp

samtools sort -T tmp -o $outBase.aligned.2.clean.shifted.bam $outBase.aligned.2.clean.shifted.bam;
samtools index $outBase.aligned.2.clean.shifted.bam;

##Splitting

###This step is done using the pipeline from Rockefeller University
###The packages are available only on R/3.4.0 so it uses a specific conda environment 
###Make sure that the packages are installed before running

conda activate r_chipqc_env;
cd $bamcleanDir
R CMD BATCH /data/rama/labMembers/gg699/atac_pipeline_test/march/Rcode_atac_2.r;
conda deactivate;


##BigWigs with deeptools

conda activate atac_chip_pipeline;
bamCoverage -b $bamsplitDir/$outBase.diNuc.bam --normalizeUsing CPM --outFileFormat bigwig --binSize 1 -o $bwigDir/$outBase.diNuc.CPM.bw
bamCoverage -b $bamsplitDir/$outBase.monoNuc.bam --normalizeUsing CPM --outFileFormat bigwig --binSize 1 -o $bwigDir/$outBase.monoNuc.CPM.bw
bamCoverage -b $bamsplitDir/$outBase.nucFree.bam --normalizeUsing CPM --outFileFormat bigwig --binSize 1 -o $bwigDir/$outBase.nucFree.CPM.bw
conda deactivate;

###Prepare bedfile for HOMER

cd $peakDir
sed 's/"//g' bed_promoters.txt > bed_promoters.bed
dos2unix bed_promoters.bed
tr [:blank:] \\t < bed_promoters.bed >  promoters.bed

bedtools intersect -a $outBase.shifted_peaks.bed -b promoters.bed > $outBase.shifted_peaks_promoters.bed -wa



##Last QC plots and footprinting
###Done in my R version 3.5.1
###You need to provide a motif of interest or a list of motifs modifying the Rcode_atac_3.r file

cd $bamcleanDir
conda activate r_3_6_env;
R CMD BATCH /data/rama/labMembers/gg699/atac_pipeline_test/march/Rcode_atac_3.r;
conda deactivate;



##Cleaning intermediate files 

rm $bamcleanDir/$outBase.filtered.3.bam
rm $bamcleanDir/$outBase.marked.2.bam
rm $bamcleanDir/$outBase.sorted.1b.bam
rm $bwaDir/$outBase.aligned.1.bam










