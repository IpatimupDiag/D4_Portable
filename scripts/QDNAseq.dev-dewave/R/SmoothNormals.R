#########################################################################/**
# @RdocFunction SmoothNormals
#
# @alias smoothobject,NoWaves-method
#
# @title "Perform a loess fitting on normal copynumber data"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{A @see "data.frame" object with log2 \code{copynumber}
#         data from normal DNA.}
# }
#
# \value{
#     Returns a @see dataframe of loess smoothed normal log2 copynumber data, 
#     to be used as input for @see "CorrectTumors".
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
setMethod('SmoothNormals', definition=function(object,bandwidth=1){
      CGHNormal <- object    
  
      fit.loess <- function(object,iterations=4,window=50,nchr=22){
        #object <- MA;iterations=4;window=50;nchr=22
        nsample <- ncol(object$M)
        nchroms <- nchr   
        object$fitted.loess <- matrix(NA,ncol=nsample,nrow=nrow(object$M[object$genes$Chr <= nchroms,]))  
        # A step to augment the data with an additional 50 points at the beginning and end of each chromosome    
        y.pre.sam <- matrix(NA,ncol=nsample,nrow=nrow(object$M[object$genes$Chr <= nchroms,])+(nchroms*100))
        weights.pre.sam <- matrix(NA,ncol=nsample,nrow=nrow(object$M[object$genes$Chr <= nchroms,])+(nchroms*100))                       
        for (j in 1:nsample) {
          y.pre.inter <- vector()
          for (k in 1:nchroms){
            y.pre <- object$M[object$genes$Chr==k, j]
            
            MAD <- mad(y.pre,na.rm=T)/4
            begin.rand <- rnorm(50,median(y.pre,na.rm=T),MAD)
            end.rand <- rnorm(50,median(y.pre,na.rm=T),MAD)
            
            new.y.pre <- c(begin.rand,y.pre,end.rand)
            y.pre.inter <- c(y.pre.inter,new.y.pre)
          }
          y.pre.sam[,j] <- y.pre.inter
        }
        
        for (j in 1:nsample) {
          weights.pre.inter <- vector()
          for (k in 1:nchroms){
            weights.pre <- object$weights[object$genes$Chr==k, j]
            
            begin.w <- rep(1,50)
            end.w <- rep(1,50)
            
            new.weights.pre <- c(begin.w,weights.pre,end.w)
            weights.pre.inter <- c(weights.pre.inter,new.weights.pre)
          }
          weights.pre.sam[,j] <- weights.pre.inter
        }
        
        x.pre <- matrix(NA,ncol=2,nrow=nrow(object$M[object$genes$Chr <= nchroms,])+(nchroms*100))
        
        x.pre.inter <- vector()
        x.pre.inter.two <- vector()
        
        for (k in 1:nchroms){
          x.pre.1 <- c(rep(k,50),object$genes$Chr[object$genes$Chr==k],rep(k,50))
          x.pre.inter <- c(x.pre.inter,x.pre.1)
          
          
          x.pre.2 <- c(seq(min(object$genes$Position[object$genes$Chr==k])-5000000,min(object$genes$Position[object$genes$Chr==k])-100000,100000),
                       object$genes$Position[object$genes$Chr == 
                                               k], seq(max(object$genes$Position[object$genes$Chr == 
                                                                                   k]) + 100000, max(object$genes$Position[object$genes$Chr == 
                                                                                                                             k]) + 5000000, 100000))
          x.pre.inter.two <- c(x.pre.inter.two,x.pre.2)
        }
        
        x.pre[,1] <- x.pre.inter
        x.pre[,2] <- x.pre.inter.two
        
        colnames(x.pre) <- c("Chr","Position")
        
        x.pre <- as.data.frame(x.pre)
        
        for (j in 1:nsample) {
          print(paste("Loess fitting for profile:",j))
          for (k in 1:nchroms) {
            y <- y.pre.sam[x.pre$Chr==k, j]
            x <- x.pre$Position[x.pre$Chr==k]
            w <- weights.pre.sam[x.pre$Chr==k, j]
            span <- window/(length(x))
            iterations <- 4
            object$M[object$genes$Chr==k, j] <- loessFit(y, x, w, span = span, 
                                                         iterations = iterations)$residuals[-(c(1:50,(length(y)-49):length(y)))]
            object$fitted.loess[object$genes$Chr[object$genes$Chr <= nchroms] ==k, j] <- loessFit(y, x, w, span = span, 
                                                                                                  iterations = iterations)$fitted[-(c(1:50,(length(y)-49):length(y)))]
          }
        }
        return(object)  #error corrected 2010-03-09
      }
      
      fitlo <- function(data,wind=20,n_ann=3){
        #data <- CGHNormal;wind=20;n_ann=3
        colnamesdata <- colnames(data)
        object <- list()
        object$M <- data[,-(1:n_ann)] ### this creates a matrix, M, which contains all of the log2 ratios
        object$genes <- matrix(NA,ncol=5,nrow=nrow(data)) ## this object will contain information about the clones
        object$genes <- as.data.frame(object$genes)
        object$genes[,1] <- data[,1]
        object$genes[,2] <- data[,2]
        object$genes[,3] <- data[,3]
        object$genes[,4] <- data[,3]
        object$genes[,5] <- data[,3] 
        object <- new("MAList",object) ## putting the object into the form of an MA list
        colnames(object$M) <- names(data)[-(1:n_ann)]
        colnames(object$genes) <- c("Clone","Chr","Start","End","Position")
        MA <- object
        MA$design <- rep(1,ncol(MA$M))
        ord <- order(MA$genes$Chr, MA$genes$Start)
        MA <- MA[ord,]
        MA$weights <- matrix(NA,nrow=nrow(MA$M),ncol=ncol(MA$M))
        for (i in 1:ncol(MA$M)){
          MA$weights[,i] <- ifelse(abs(MA$M[,i]) > 0.3, 0, 1)
        }
        set.seed(1) # this sets the seed equal to 1
        MA.L <- fit.loess(MA,iterations=4,window=wind,nchr=22)
        colnames(MA.L$fitted.loess) <- colnames(MA.L$M)
        MA.fit <- MA.L$fitted.loess
        return(MA.fit)
      }
      
      chromos <- CGHNormal[,2]
      wh_no_24 <- which(chromos<=22)
      bw <- max(20,round(20*bandwidth*nrow(CGHNormal)/44000))
      CGHNormal <- CGHNormal[wh_no_24,]
      if(mean(CGHNormal[,4],na.rm = TRUE) > 100000) {n_end <- 4} else {n_end <- 3}
      
      if(length(which(is.na(CGHNormal)))>0){
        print("Missings detected. Started imputation")
        data_imp <- impute.knn(as.matrix(CGHNormal[,-(1:n_end)]))
        if(mode(data_imp)=="list") CGHNormal <- data.frame(CGHNormal[,1:n_end],data_imp[[1]]) else CGHNormal <- data.frame(CGHNormal[,1:n_end],data_imp)
      } else {print("No missings detected")}
      
      NormalsSm0 <- fitlo(data=CGHNormal,wind=bw,n_ann=n_end)
      NormalsSm0Ann <- data.frame(CGHNormal[,1:n_end],NormalsSm0)
      return(NormalsSm0Ann)
    
})