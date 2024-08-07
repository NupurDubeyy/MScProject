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
#download.file(url="https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz", destfile="Homo_sapiens.GRCh38.112.chr.gtf.gz", timeout=300)
#GTFfile = "/data4/msc23104469/dexseq/Homo_sapiens.GRCh38.112.chr.gtf"
#txdb = makeTxDbFromGFF("GTFfile", format = "gtf")

#file_to_be_removed = '/data4/msc23104469/dexseq/dexseq_breast/Homo_sapiens.GRCh38.112.chr.gtf.gz'
#file.remove(file_to_be_removed)
#print("The file 'Homo_sapiens.GRCh38.112.chr.gtf.gz'is removed")

# Define the URL and the destination file path
url <- "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
destfile <- "Homo_sapiens.GRCh38.112.gtf.gz"

# Use the system function to call wget
system(paste("wget", url, "-O", destfile))

# Load the GenomicFeatures package
library(GenomicFeatures)

# Specify the full path to the downloaded GTF file
gtf_file <- file.path(getwd(), destfile)

# Check if the file exists
if (file.exists(gtf_file)) {
  # Create a TxDb object from the GTF file
  txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")
} else {
  stop("The specified GTF file does not exist.")
}

flattenedAnnotation = exonicParts( txdb, linked.to.single.gene.only=TRUE )
names(flattenedAnnotation) = sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)
print("flattenedAnnotation:")
head(names(flattenedAnnotation))
## ----Counting reads------------------------------------------------------------

print("Counting reads...")
print("reading bam folder...")
bam_folder <- "/data4/msc23104469/alignment/alignment_BC/bams"

bam_files <- list.files(bam_folder, pattern = "\\.bam$", full.names = TRUE)
print("bam files list complete:")
bam_files <- c(bam_files)
head(bam_files)

bamFilesobject = BamFileList( bam_files, asMates = TRUE )
# add this yieldSize = 2e5 if script cancels, between the 2 args
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

print("running estimateDispersions...")
dxd = estimateDispersions( dxd, formula = formulaFullModel, BPPARAM = BPPARAM )
saveRDS(dxd, file = "dxd_estimateDispersions.Rds")

print("running testForDEU...")
dxd = testForDEU( dxd, reducedModel = formulaReducedModel_PE, fullModel = formulaFullModel, BPPARAM = BPPARAM )
saveRDS(dxd, file = "dxd_testForDEU.Rds")

print("running estimateExonFoldChanges...")
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM = BPPARAM)
saveRDS(dxd, file = "dxd_estimateExonFoldChanges.Rds")

print("extracting DEXSeqResults...")
res = DEXSeqResults( dxd )
saveRDS(res, "dxr_v4.rds")

print("Summarizing results on gene level...")
pgq <- perGeneQValue(res, p = "pvalue")

print("dxr object created. Writing html report")
DEXSeqHTML( res, FDR=0.05, color=c("#FF000080", "#0000FF80") )

