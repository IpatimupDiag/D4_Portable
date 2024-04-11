#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq analysis segmentBins
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

source("scripts/plotQDNAseq.R", echo  = !msg)

dewaved<- snakemake@input[["dewaved"]]
bin <- as.integer(snakemake@wildcards[["binSize"]])
segmented <- snakemake@output[["segmented"]]
profiles <- snakemake@params[["profiles"]]
failed<- snakemake@params[["failed"]]
min_used_reads<-snakemake@params[["minimal_used_reads"]]

copynumbersbed<-snakemake@params[["copynumbersbed"]]
segmentsbed<-snakemake@params[["segmentsbed"]]
copynumbers<-snakemake@output[["copynumbers"]]
segments<-snakemake@output[["segments"]]
bedfolder <- snakemake@params[["bedfolder"]]

log<-snakemake@log[[1]]
log<-file(log, open="wt")
sink(log, append=TRUE , split=FALSE)

# Adjust segmentation settings based on binsize
if (bin==15) {SDundo=0.75; alph=1e-15}
if (bin==30) {SDundo=0.75; alph=1e-15}
if (bin==100) {SDundo=0.10; alph=1e-20} # default = 1e-20 # 0.01 used for PELLL_FS8_a0.01
if (bin==1000) {SDundo=0.10; alph=1e-20}
# TODO: mogelijk hogere alpha nodig, om meer segmenten te krijgen: 0.01 (1e-2). Deze setting used in tmp_100kbp

# load data
QCN.fcnsd <- readRDS(dewaved)

QCN.fcnsds <- segmentBins(QCN.fcnsd[,QCN.fcnsd$used.reads > min_used_reads ], undo.splits='sdundo', undo.SD=SDundo, alpha=alph, transformFun="sqrt") # gives 'Performing segmentation: NA
QCN.fcnsdsn <- normalizeSegmentedBins(QCN.fcnsds)
saveRDS(QCN.fcnsdsn, segmented)

##############################################################################################################
# Create profiles
##############################################################################################################

plotQDNAseq(QCN.fcnsdsn, profiles)

#create output for failed samples - for snakemake compatibility.
littledata<-colnames(QCN.fcnsd[,QCN.fcnsd$used.reads <= min_used_reads ])
if(length(littledata>0)){for(file in littledata){file.create(paste(profiles, file,".png",sep=""))
file.create(paste(bedfolder, file,"-copynumbers.bed",sep=""))
file.create(paste(bedfolder, file,"-segments.bed",sep=""))
}}

write.table(littledata, file=failed )

##############################################################################################################
# Create IGV objects and bedfiles from readcounts
##############################################################################################################

exportBins(QCN.fcnsdsn, copynumbers, format="igv", type="copynumber")
exportBins(QCN.fcnsdsn, segments, format="igv", type="segments")
exportBins(QCN.fcnsdsn, file=copynumbersbed, format="bed", logTransform=TRUE, type="copynumber")
exportBins(QCN.fcnsdsn, file=segmentsbed, format="bed", type="segments")
