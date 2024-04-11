#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis to normalize bins and dewave
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville)
# date: December 2017
# Changed to work in snakemake pipeline by Tjitske Los
##############################################################################################################
msg <- snakemake@params[["suppressMessages"]]
if (msg){
suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))
} else{
library(QDNAseq)
library(Biobase)
}

source("scripts/functions.R", echo = !msg)
source("scripts/plotQDNAseq.R", echo = !msg)

binReadCounts <- snakemake@input[["binReadCounts"]]
bin <- as.integer(snakemake@wildcards[["binSize"]])
corrected <- snakemake@output[["corrected"]]
profiles <- snakemake@params[["profiles"]]
chrom_filter <- c(snakemake@params[["chrom_filter"]])

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T , split=FALSE)
##############################################################################################################
# Correct, Normalize & Dewave raw data
##############################################################################################################
QRC <- readRDS(binReadCounts)

#did these steps manualy as is described in second part creationofblacklist.R for samples with chrX.

QRC.f <- applyFilters(QRC, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE , chromosomes=chrom_filter)

#featureData(QRC.f)$use[c(22232:22237)]<-F  # segment 14:22400001-22500000 t/m 14:22900001-23000000 on chr14
QRC.f <- estimateCorrection(QRC.f)
#QRC.f <- applyFilters(QRC, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE)  - used for bia-ALCL (2019-07-09 not)
#QRC.f <- estimateCorrection(QRC.f, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE) - used for bia-ALCL(2019-07-09 not)
#QRC.f <- applyFilters(QRC.f, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE, chromosome="Y") - used for bia-ALCL(2019-07-09 not)

QCN.fc <- correctBins(QRC.f)
QCN.fcn <- normalizeBins(QCN.fc)
QCN.fcns <- smoothOutlierBins(QCN.fcn)

saveRDS(QCN.fcns, corrected)
plotQDNAseq(QCN.fcns, profiles)
