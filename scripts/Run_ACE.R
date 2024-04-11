#!/usr/bin/env Rscript
##############################################################################

#Author: Tjitske logs
#date: Dec 2017

#This script is a small wrapper around ACE to work in the snakemake pipeline
##############################################################################
msg <- snakemake@params[["suppressMessages"]]
if (msg){
suppressMessages(library(QDNAseq))
suppressMessages(library(stringr))
} else{
library(QDNAseq)
library(stringr)
}

source('scripts/ACE.R', echo = !msg)

ploidies<-as.integer(snakemake@wildcards[["ploidy"]])
inputfile <-snakemake@input[["segmented"]]
outputdir<-snakemake@params[["outputdir"]]
failed <- snakemake@params[["failed"]]
log<-snakemake@log[[1]]
fitpicker <- snakemake@output[["fitpicker"]] # Insert by Erik 2021-01-27

imagetype <- snakemake@config[["ACE"]][["imagetype"]]
method<-snakemake@config[["ACE"]][["method"]]
penalty<-as.numeric(snakemake@config[["ACE"]][["penalty"]])
cap<-as.integer(snakemake@config[["ACE"]][["cap"]])
trncname<- snakemake@config[["ACE"]][["trncname"]]
printsummaries<- snakemake@config[["ACE"]][["printsummaries"]]

copyNumbersSegmented <- readRDS(inputfile)

parameters <- data.frame(options = c("ploidies","imagetype","method","penalty","cap","trncname","printsummaries"),
                         values = c(paste0(ploidies,collapse=", "),imagetype,method,penalty,cap,trncname,printsummaries))

#write.table(parameters, file=log, quote = FALSE, sep = "\t", na = "", row.names = FALSE)
write.table(parameters, file=log, quote = FALSE, sep = "\t", na = "", col.names = NA)

ploidyplotloop(copyNumbersSegmented ,outputdir , ploidies,imagetype,method,penalty,cap,trncname,printsummaries) # 1 to 3 dir.create warnings and 4 to 15  removed rows warnings, output gives NULL

#create output for failed samples - for snakemake compatibility.
failed_samples<-read.table(failed, stringsAsFactors=FALSE, header=TRUE)
if(length(failed_samples[,1]>0)){
  for(file in failed_samples[,1]){file.create(paste(outputdir, ploidies,"N/", file,"/summary_",file,".",imagetype,sep=""))}
}

