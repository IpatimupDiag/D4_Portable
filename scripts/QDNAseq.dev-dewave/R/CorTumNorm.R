#########################################################################/**
# @RdocFunction CorTumNorm
#
# @alias smoothobject,NoWaves-method
#
# @title "Creating a correlation matrix between tumor and normal copynumber data"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{A @see "data.frame" object with log2 \code{copynumber}
#         data.}
#     \item{user.calibration}{A @see "data.frame" object with loess fitted normal log2 \code{copynumber} data}
# }
#
# \value{
#     Returns a @see dataframe of with a correlation between tumor and normal copynumber data. 
# }
#
# \section{How to cite this package}{
#  Whenever using this package, please cite:
#   van de Wiel MA, Brosens R, Eilers PH, Kumps C, Meijer GA, Menten B, Sistermans E, Speleman F, Timmerman ME, Ylstra B.
#   "Smoothing waves in array CGH tumor profiles."
#   Bioinformatics. 2009 May 1;25(9):1099-104. doi: 10.1093/bioinformatics/btp132. Epub 2009 Mar 10.
# }
#
# @author "MC"
# \seealso{
#     Internally, \code{NoWaves} of the NoWaves package.
# }
#
# @keyword manip
#*/#########################################################################
setMethod('CorTumNorm', signature=c(object="data.frame",user.calibration="data.frame"),
          definition= function(object,user.calibration){
        
        CGHTumor <- object
        CGHNormal <- user.calibration
  
        if(mean(CGHNormal[,4],na.rm = TRUE) > 100000) {n_end <- 4} else {n_end <- 3}
        if(mean(CGHTumor[,4],na.rm = TRUE) > 100000) {t_end <- 4} else {t_end <- 3}
        mergeddata <- merge(CGHTumor,CGHNormal[,-c(2:n_end)],by=1,sort=F) #merges on probe name
        nnorm <- ncol(CGHNormal)
        ntum <- ncol(CGHTumor)
        CGHNormal <- mergeddata[, c(1:t_end, (ntum + 1):(nnorm + ntum - n_end))]
        CGHTumor <- mergeddata[,1:ntum,drop=FALSE]
        cormat <- signif(cor(cbind(CGHTumor[,-(1:t_end)],CGHNormal[,-(1:n_end)]),use="pairwise.complete.obs"),2)
        cormat_part <- signif(cormat[1:(ntum-t_end),-(1:(ntum-t_end))],2)
        return(cormat_part) #displays correlation matrix
        
    })