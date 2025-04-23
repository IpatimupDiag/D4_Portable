#!/usr/bin/env Rscript

#TO BE Adapted TO SANKEMAKE:
##############################################################################################################
# script for Clonality
# Adapted from /ccagc/lib/pipelines/QDNAseq/QDNAseq.R (Daoud Sie) -and MM_QDNAseq (Matias Mendeville) # TO UPDATE !!!
# date: March 2024
# Changed to work in snakemake pipeline by José Pedro Matos
##############################################################################################################

msg <- snakemake@params[["suppressMessages"]]
.libPaths(c("../.R/lib/", .libPaths()))
#print(.libPaths())


#print('////////////////////////////')

library(devtools)

if(!require("QDNAseq.dev"))
	devtools::install_github("tgac-vumc/QDNAseq.dev", ref="clonality",force=TRUE)

library(QDNAseq.dev)
library(Biobase)
library(DNAcopy)
library(gridExtra)
library(grid)
library(Clonality)


#source("scripts/d4-pipeline-master/bin/clonality.R", echo  = !msg) # TO UPDATE !!! # WHAT IT DOES ???
#============================================================================
pjct <- snakemake@params[["projname"]]
RDS <- snakemake@input[["RDS"]]
ref <- snakemake@params[["reference"]]
profiles_dir <- snakemake@output[["profiles"]]
report_png <- snakemake@output[["report_png"]]
report_bmp <- snakemake@output[["report_bmp"]]
#============================================================================

#print(RDS) #DEBUG
#print(',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,')

GoldStandard_Data <- read.delim('scripts/d4-pipeline-master/data/TRACERx_GoldStandard_Clonality_Statistics.txt')
GMM_model <- readRDS('scripts/d4-pipeline-master/data/GMM_model.Rds')
refData <- read.table(ref,skip=0,fill=TRUE, header=TRUE, sep="\t",quote="")


duplicated(refData[,3]) -> exc
CNA(refData[!exc,c(-1,-2,-3,-4)], refData[!exc,2], refData[!exc,4] + refData[!exc,3] / 2) -> refCNA

data <- readRDS(RDS)

sub("(_S\\d+)?(_L\\d{3})?_R\\d{1}_\\d{3}(\\.f(ast)?q\\.gz)?$", "", sampleNames(data)) -> samples

if (file.exists("key")) {
	key <- read.delim("key", sep="\t", header=F, stringsAsFactors=F)
	matched <- match(samples, key[,1])
	data <- data[,!is.na(matched)]
	samples <- key[matched[!is.na(matched)],2]
}

sub("CNP\\d{2}-\\d{3}-", "", samples) -> samples
sub(paste(pjct, "-", sep=""), "", samples) -> samples
sampleNames(data) <- samples
#pData(data)$name <- samples

if (sum(duplicated(samples)) > 0) {
	QRC <- poolRuns(data, samples=samples, force=TRUE)
	QCN.fcnss <- segmentBins(QRC)
	QCN.fcnssn <- normalizeSegmentedBins(QCN.fcnss)
	data <- QCN.fcnssn
}
CNA(assayDataElement(data, "copynumber"), chromosomes(data), bpstart(data)) -> dataCNA

clonality.analysis(dataCNA, ptlist=rep(pjct, ncol(data)), refdata=refCNA) -> llrData

clonalityTest(data, sbjctLst = rep(pjct, ncol(data)), refdata=refCNA, llrData=llrData) -> cln

clonTab <- cln$clonTab
cor <- 1 - clonTab$cor
clonTab$cor <- cor

# Determine GMM clonality 
model_prediction <- predict(GMM_model,newdata = clonTab, type = 'response')
Non_clonal_prob <- model_prediction
Clonal_prob <-  1 - Non_clonal_prob
# Check if sample is potential outlier: GMM and twometric don't match
if((clonTab$cor > 0.54 & clonTab$llr2 > 0) & Non_clonal_prob > 0.5){
    outlier <- T
    main <- 'This sample is an outlier and possibly non-clonal'
}else if((clonTab$llr2 < -5 | clonTab$cor < 0.45) & Non_clonal_prob < 0.5){
    outlier <- T
    main <- 'This sample is an outlier and possibly clonal'
}else{
    outlier <- F
    main = ''
}


#Save profiles and reports:
png(report_png, height=297, width=210, unit="mm", res=150) #png(paste(pjct, ".report.png", sep="") , height=297, width=210, unit="mm", res=150)
layout(1:2)
clonalityReport(cln, labels=cln$clonTab$combination.named)
# Add gold standard points in the background
points(GoldStandard_Data$llr2,GoldStandard_Data$cor,  pch=ifelse(as.integer(as.factor(GoldStandard_Data$clonality)) == 1,4,6), cex=1, col =rgb(0, 0, 0, 0.15))
# add outlier text if any
text(x= 50, y = -0.8,main, col = 'red')
grid.table(clonTab[,c(1,3,4,5)], vp=viewport(x=unit(0.5, "npc"), y=unit(0.25, "npc")))
dev.off()

bmp(report_bmp, height=297, width=210, unit="mm", res=150)
layout(1:2)
clonalityReport(cln, labels=cln$clonTab$combination.named)
# Add gold standard points in the background
points(GoldStandard_Data$llr2,GoldStandard_Data$cor,  pch=ifelse(as.integer(as.factor(GoldStandard_Data$clonality)) == 1,4,6), cex=1, col =rgb(0, 0, 0, 0.15))
# add outlier text if any
text(x= 50, y = -0.8,main, col = 'red')
grid.table(clonTab[,c(1,3,4,5)], vp=viewport(x=unit(0.5, "npc"), y=unit(0.25, "npc")))
dev.off()

png(paste(profiles_dir,paste(pjct, ".profiles.%02d.png", sep=""), sep=""), height=297, width=210, unit="mm", res=150)
layout(1:3)
plot(data)
dev.off()

bmp(paste(profiles_dir,paste(pjct, ".profiles.%02d.bmp", sep=""), sep=""), height=297, width=210, unit="mm", res=150)
layout(1:3)
plot(data)
dev.off()

for (i in 1:ncol(data)) {
	sn <- sampleNames(data)[i]
	png(paste(profiles_dir,paste(pjct, sn, "png", sep="."), sep=""), width=960, height=480)
	plot(data[,i])
	dev.off()	
	bmp(paste(profiles_dir,paste(pjct, sn, "bmp", sep="."), sep=""), width=960, height=480)
	plot(data[,i])
	dev.off()
}

