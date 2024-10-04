#

# Create new calibration profiles for deWave from normals.
# All variables with '!!!'  as coment, chanhe the bins size value.

#
library(gtools)
library(impute)
library(MASS)
library(limma)
library(QDNAseq)
library(Biobase)

sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm,":")           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
}
sourceDir("R/") # get deWave local lib

#Normals input directory

## ionTurrent
#in_dir <- "/home/jpparrachadematos/new_D4_dir/normalsData_deWave/IPATIMUP_normals/IPATIMUP_normals_packALL/"

## illumina
#in_dir <- "/home/jpparrachadematos/new_D4_dir/normalsData_deWave/TGAC_normals/TGAC_normals_pck/"

## mix of both
in_dir <- "/home/jpparrachadematos/new_D4_dir/normalsData_deWave/mix_normals/"

#Bin size !!!
# Change for the disired bin size.
bin_size <- 1000 #!!!

#Calibration profiles directory
out_dir <- "data_mix/" #!!!

#start

bins <- getBinAnnotations(binSize=bin_size)

print("Get the normals:")

readCounts <- binReadCounts(bins, path=in_dir)

print("Filtering the normals:")

readCountsFiltered <- applyFilters(readCounts)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)

print(copyNumbers)

print("Making the calibration profiles:")

# NormalCalibrationSet_100kb !!!
# Change the bins size in the name of the variable, this is due to the name of the dataset beeing the same as the variable.
NormalCalibrationSet_1000kb <- createNormalCalibration(copyNumbers) #!!!

#save calibration profiles

calibration_set_name <- paste("NormalCalibrationSet_", bin_size, "kb.rda",sep="")
calibration_dir <- paste(out_dir, calibration_set_name,sep="")

# NormalCalibrationSet_<n>kb !!!
# Also change were.
save(NormalCalibrationSet_1000kb,file=calibration_dir) #!!!

print(paste("Done!", calibration_set_name))
