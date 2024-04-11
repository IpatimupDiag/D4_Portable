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

# for dewaving:
suppressMessages(library(limma))
suppressMessages(library(gtools))
suppressMessages(library(impute))
suppressMessages(library(MASS))
} else{
library(QDNAseq)
library(Biobase)

# for dewaving:
library(limma)
library(gtools)
library(impute)
library(MASS)
}


source("scripts/functions.R", echo = !msg)
source("scripts/plotQDNAseq.R", echo = !msg)
sourceDir(snakemake@config[["QDNAseq"]][["dewave_dir"]], echo = !msg)
load(snakemake@params[["dewave_data"]], verbose = !msg)

bin <- as.integer(snakemake@wildcards[["binSize"]])
corrected <- snakemake@input[["corrected"]]
profiles <- snakemake@params[["profiles"]]
dewaved <- snakemake@output[["dewaved"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=T , split=FALSE)
##############################################################################################################
# Correct, Normalize & Dewave raw data
##############################################################################################################

QCN.fcns <- readRDS(corrected)

if(bin!=1000){
QCN.fcns <- dewaveBins(QCN.fcns)
}

saveRDS(QCN.fcns, dewaved)
plotQDNAseq(QCN.fcns, profiles)
