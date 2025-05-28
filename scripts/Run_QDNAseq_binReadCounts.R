#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis from fastq to called readcounts
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville)
# date: December 2017
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
msg <- snakemake@params[["suppressMessages"]]
if (msg){
suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))
suppressMessages(library(R.cache))
} else{
library(QDNAseq)
library(Biobase)
library(R.cache)
}


setCacheRootPath(path="../.Rcache")

bam <- snakemake@input[["bams"]]
bin <- as.integer(snakemake@wildcards[["binSize"]])
genome <- snakemake@params[["genome"]]
binReadCounts <- snakemake@output[["binReadCounts"]]

selected_bin_100K <- snakemake@params[["CUSTOM_BINS_100K"]]
selected_bin_1000K <- snakemake@params[["CUSTOM_BINS_1000K"]]

##############################################################################################################
# Get bin annotations and bin read counts
##############################################################################################################
#custom bin annotations for chrX used - did this whole script by hand for 100kbp See creationofblacklist.R


###
if (genome == "hg19") {
        print("Using Custom bins")
        print("")
        if (bin == 100){
                print("Bin_size: 100")
                print("")
                bins <- readRDS(selected_bin_100K)
        } else if (bin == 1000) {
                print("Bin_size: 1000")
                print("")
                bins <- readRDS(selected_bin_1000K)
        }else{
                print("Using the original QDNAseq bins")
                print("The custom bins are limited to 100Kbp and 1000Kbp")
                print("")
                bins <- getBinAnnotations(bin, genome=genome)
     }	
}else{
        print("Using the original QDNAseq bins")
        print("The custom bins are limited to hg19 reference")
        print("")
        bins <- getBinAnnotations(bin, genome=genome)
}

###

QRC <- binReadCounts(bins, bamfiles=bam, cache=TRUE)

# Since it will run from .bam files:
#sub("(_[ACGT]+)?(_S\\d+)?(_L\\d{3})?_R\\d{1}_\\d{3}(\\.f(ast)?q\\.gz)?$", "", sampleNames(QRC)) -> samples
sub("(_[ACGT]+)?(_S\\d+)?(_L\\d{3})?_R\\d{1}_\\d{3}(\\.f(ast)?q\\.gz)?$", "", sampleNames(QRC)) -> samples
#sampleNames(QRC)) -> samples #TO BE TESTED !!!

if (sum(duplicated(samples)) > 0) {
        QRC <- poolRuns(QRC, samples=samples, force=TRUE)
}
saveRDS(QRC, binReadCounts)
