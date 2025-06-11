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

custom_hg19_100Kbp <- snakemake@params[["CUSTOM_BINS_100K"]]
custom_hg19_1000Kbp <- snakemake@params[["CUSTOM_BINS_1000K"]]

##############################################################################################################
# Get bin annotations and bin read counts
##############################################################################################################
#custom bin annotations for chrX used - did this whole script by hand for 100kbp See creationofblacklist.R

#bins <- getBinAnnotations(bin, genome=genome)

#Use of custom hg19 bin, size 100Kbp and 1000Kbp:
if (genome=="hg19") {
        if (bin==100) {
                print("Using custom 100Kpb hg19 bins !")
                bins <- readRDS(custom_hg19_100Kbp)
        } else if (bin==1000){
                print("Using custom 1000Kpb hg19 bins !")
	        bins <- readRDS(custom_hg19_1000Kbp)
	}
} else {
        print("Using original bins !")
        bins <- getBinAnnotations(bin, genome=genome)
}




QRC <- binReadCounts(bins, bamfiles=bam, cache=TRUE)


# Since it will run from .bam files:
#sub("(_[ACGT]+)?(_S\\d+)?(_L\\d{3})?_R\\d{1}_\\d{3}(\\.f(ast)?q\\.gz)?$", "", sampleNames(QRC)) -> samples
sub("(_[ACGT]+)?(_S\\d+)?(_L\\d{3})?_R\\d{1}_\\d{3}(\\.f(ast)?q\\.gz)?$", "", sampleNames(QRC)) -> samples
#sampleNames(QRC)) -> samples #TO BE TESTED !!!

if (sum(duplicated(samples)) > 0) {
        QRC <- poolRuns(QRC, samples=samples, force=TRUE)
}
saveRDS(QRC, binReadCounts)
