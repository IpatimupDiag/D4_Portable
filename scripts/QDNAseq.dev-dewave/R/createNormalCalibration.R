#########################################################################/**
# @RdocFunction createNormalCalibration
#
# @alias smoothobject,QDNAseqCopyNumbers-method
#
# @title "Creating a calibration set for dewaving by using a set of copynumber data from normal DNA"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{A @see "QDNAseqCopyNumbers" object with \code{copynumber}
#         data.}
#     \item{force}{Running this function will remove possible segmentation and
#         calling results. When they are present, running requires specifying
#         \code{force} is @TRUE.}
# }
#
# \value{
#     Returns a @see dataframe of loess smoothed normal copynumber profiles, 
#     to be used as calibration set for the dewaving of tumor copynumber profiles
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# copyNumbers <- correctBins(readCountsFiltered)
# NormalCalibrationSet <- createNormalCalibration(setofnormals)
# copyNumbersDewaved <- dewaveBins(copyNumbers, calibration=NormalCalibrationSet)
# }
#
# @author "MC"
#
# @keyword manip
#*/#########################################################################
setMethod('createNormalCalibration', signature=c(object='QDNAseqCopyNumbers'),
          definition=function(object,force=FALSE) {
            
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Validate arguments
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Argument "object":
            if (!force && ("segmented" %in% assayDataElementNames(object)))
              stop("Data has already been segmented. Re-normalizing will ",
                   "remove segmentation (and possible calling) results. ",
                   "Please specify force=TRUE, if you want this.")
            
            if ("segmented" %in% assayDataElementNames(object))
              assayDataElement(object, "segmented") <- NULL
            if ("calls" %in% assayDataElementNames(object)) {
              assayDataElement(object, "calls") <- NULL
              assayDataElement(object, "probloss") <- NULL
              assayDataElement(object, "probnorm") <- NULL
              assayDataElement(object, "probgain") <- NULL
              if ("probdloss" %in% assayDataElementNames(object))
                assayDataElement(object, "probdloss") <- NULL
              if ("probamp" %in% assayDataElementNames(object))
                assayDataElement(object, "probamp") <- NULL
            }


            #Create calibration dataframe
            data <- assayData(object)
            feature <- featureData(object)
            chromosome <- feature$chromosome
            chromosome <- gsub("X", 23, chromosome)
            chromosome <- gsub("Y", 24, chromosome)
            chromosome <- as.numeric(chromosome)
            
            start <- feature$start
            end <- feature$end
            
            copynumber <- assayData(object)
            copynumber <- copynumber$copynumber
            copynumber <- copynumber[,-1]
            copynumber[copynumber==0] <- 1
            copynumber <- log2(copynumber)
            
            dataset <- data.frame(Name=1:nrow(copynumber), Chromosome=as.vector(chromosome), Start=start, End=end)
            dataset <- cbind(dataset, copynumber)
            
            smooth_input <- na.omit(dataset)
            
            object <- SmoothNormals(smooth_input,bandwidth=1)
            
            object
            
          })

# EOF