#!/usr/bin/env Rscript
##############################################################################

#Author: Tjitske logs
#date: Feb 2018

#This script is a small wrapper around ACE to perform the postanalysisloop in the snakemake pipeline
##############################################################################
msg <- snakemake@params[["suppressMessages"]]
if (msg){
suppressMessages(library(QDNAseq))
} else{
library(QDNAseq)
}

source('scripts/ACE.R', echo = !msg)

## Erik
#copyNumbersSegmented <- snakemake@output[["ACE_post"]]
##

inputfile <-snakemake@input[["segmented"]]
outputdir<-snakemake@params[["outputdir"]]
failed <- snakemake@params[["failed"]]
fitpickertable<-snakemake@input[["fitpicker"]]
ploidies<-as.integer(snakemake@wildcards[["ploidy"]])

imagetype <- snakemake@config[["ACE"]][["imagetype"]]
trncname<- snakemake@config[["ACE"]][["trncname"]]

copyNumbersSegmented <- readRDS(inputfile)

postanalysisloop(copyNumbersSegmented, modelsfile=fitpickertable, imagetype=imagetype, outputdir=outputdir, trncname=trncname)

# if Rplots.pdf exists rename from QDNAseq-snakemake directory to outputdir
if (file.exists("./Rplots.pdf")){
    file.rename(from = "./Rplots.pdf", to = paste0(outputdir,"all_samples.pdf"))
}

failed_samples<-read.table(failed, stringsAsFactors=FALSE, header=TRUE)
if(length(failed_samples[,1]>0)){
  for(file in failed_samples[,1]){file.create(paste(outputdir,"segmentfiles/",file,"_segments.tsv",sep=""))}
}
