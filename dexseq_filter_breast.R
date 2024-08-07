#!/usr/bin/Rscript

## ----warning=FALSE------------------------------------------------------------
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(AnnotationDbi)
library(Biobase)
library(DESeq2)
library(RColorBrewer)
library(Biostrings)
library(XVector)
library(GenomicFeatures)
library(DEXSeq)
library(SummarizedExperiment)
library(MatrixGenerics)
library(matrixStats)
library(GenomicAlignments)
library(BiocParallel)
library(S4Vectors)
library(GenomeInfoDb)
library(BiocGenerics)
library(Rsamtools)

## ----Preparing the annotation-------------------------------------------------------------------------

print("Preparing the annotation...")
#download.file("https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz", destfile="Homo_sapiens.GRCh38.112.chr.gtf.gz")

txdb = makeTxDbFromGFF(file="/data4/msc23104469/dexseq/Homo_sapiens.GRCh38.112.chr.gtf", format="gtf")

#file_to_be_removed = '/data4/msc23104469/dexseq/dexseq_filter_l/Homo_sapiens.GRCh38.112.chr.gtf.gz'
#file.remove(file_to_be_removed)
#print("The file 'Homo_sapiens.GRCh38.112.chr.gtf.gz' removed")

flattenedAnnotation = exonicParts( txdb, linked.to.single.gene.only=TRUE )
names(flattenedAnnotation) = sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)
print("flattenedAnnotation:")
head(names(flattenedAnnotation))
## ----Counting reads------------------------------------------------------------

print("Counting reads...")
print("reading bam folder...")
bam_folder <- "/data4/msc23104469/alignment/alignment_BC/bams"

bam_files <- list.files(bam_folder, pattern = "\\.bam$", full.names = TRUE)
print("bam files list completed:")
bam_files <- c(bam_files)
head(bam_files)

bamFilesobject = BamFileList( bam_files, asMates = TRUE )
# add this yieldSize = 2e5 if script cancels, between the 2 args
print(bamFilesobject)
print("bamFilesobject completed")

seqlevelsStyle(flattenedAnnotation) = "NCBI" 
BPPARAM <- SerialParam()

se <- NULL
se = summarizeOverlaps(flattenedAnnotation, bamFilesobject, singleEnd=FALSE, fragments=TRUE, ignore.strand=TRUE, BPPARAM = BPPARAM)
# Check if 'se' was created successfully and print "Success!" if so
if (!is.null(se)) {
  print("Success!")
} else {
  print("Error: The object 'se' was not created.")
}

saveRDS(se, "se.rds")

## -----Building a DEXSeqDataSet------------------------------------------------------------------------

print("Building a DEXSeqDataSet...")
samplesheet_path = "/data4/msc23104469/dexseq/breast_cancer_metadata.csv"
samplesheet <- read.csv(samplesheet_path)
head(samplesheet)

colData(se)$condition = factor(c(samplesheet$Condition))
colData(se)$patient = factor(c(samplesheet$Patient))

formulaFullModel    =  ~ sample + exon + patient:exon + condition:exon
formulaReducedModel_PE =  ~ sample + exon + patient:exon
formulaReducedModel_CE =  ~ sample + exon + condition:exon
dxd = DEXSeqDataSetFromSE( se, design= ~ sample + exon + patient:exon + condition:exon )

# Print the DEXSeqDataSet object
print(dxd)

# Display summary information
summary(dxd)

print("DEXSeqDataSeqFromSE finished...")
print("writing DEXSeqDataSetFromSE output...")
saveRDS(dxd, file = "dxd2.Rds")

## ----para1,cache=TRUE,results='hide', eval=FALSE------------------------------
print("running estimateSizeFactors...")
dxd = estimateSizeFactors( dxd )
saveRDS(dxd, file = "dxd_estimateSizeFactors.Rds")

print("filtering out exons with low count variance...")
norm_counts <- featureCounts(dxd, normalized = TRUE)
filter_bin_variance <- 0.029
vars <- apply(log2(norm_counts + 1), 1, var)
keep <- which(vars >= filter_bin_variance)
# subset dxd to only keep exons with variance > threshold
dxd_filtered <- dxd[keep, ]
print("writing dxd_filtered to file...")
saveRDS(dxd_filtered, file = "dxd_filtered_variance_0.029.Rds")

print("running estimateDispersions...")
dxd_filtered = estimateDispersions( dxd_filtered, formula = formulaFullModel, BPPARAM = BPPARAM )
saveRDS(dxd_filtered, file = "dxd_estimateDispersions.Rds")

print("running testForDEU...")
dxd_filtered = testForDEU( dxd_filtered, reducedModel = formulaReducedModel_PE, fullModel = formulaFullModel, BPPARAM = BPPARAM )
saveRDS(dxd_filtered, file = "dxd_testForDEU.Rds")

print("running estimateExonFoldChanges...")
dxd_filtered = estimateExonFoldChanges( dxd_filtered, fitExpToVar="condition", BPPARAM = BPPARAM)
saveRDS(dxd_filtered, file = "dxd_estimateExonFoldChanges.Rds")

print("extracting DEXSeqResults...")
res = DEXSeqResults( dxd_filtered )
saveRDS(res, "dxr_v4.rds")

print("Summarizing results on gene level...")
pgq <- perGeneQValue(res, p = "pvalue")

print("dxr object created. Writing html report")
DEXSeqHTML( res, FDR=0.05, color=c("#FF000080", "#0000FF80") )

