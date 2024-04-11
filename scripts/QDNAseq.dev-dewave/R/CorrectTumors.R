#########################################################################/**
# @RdocFunction CorrectTumors
#
# @alias smoothobject,NoWaves-method
#
# @title "Correcting tumor copynumber profiles based correlated genomic waves found in normal copynumber data"
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
#     Returns a @see dataframe of loess smoothed tumor copynumber profiles, 
#     to be used as input for a @see "QDNAseqCopyNumbers" object.
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
setMethod('CorrectTumors', signature=c(object="data.frame",user.calibration="data.frame"),
          definition= function(object,user.calibration,bandwidth=1,thr=2.5,
                               applyridgereg = TRUE,normalization="none"){
            
            CGHTumor <- object
            CGHNormalSmooth <- user.calibration


            #CGHTumor<-CGHTumor;CGHNormalSmooth<-NormalsSmooth;bandwidth=1;thr=2.5;applyridgereg = T;normalization="none"
            if(mean(CGHNormalSmooth[,4],na.rm = TRUE) > 100000) {n_end <- 4} else {n_end <- 3}
            if(mean(CGHTumor[,4],na.rm = TRUE) > 100000) {t_end <- 4} else {t_end <- 3}
            
            if(length(which(is.na(CGHTumor)))>0){
              print("Missings detected. Started imputation")
              data_imp <- impute.knn(as.matrix(CGHTumor[,-(1:t_end),drop=FALSE]))
              if(mode(data_imp)=="list") CGHTumor <- data.frame(CGHTumor[,1:t_end],data_imp[[1]]) else CGHTumor <- data.frame(CGHTumor[,1:t_end],data_imp)
            } #else {print("No missings detected")}
            
            bw <- max(20,round(20*bandwidth*nrow(CGHNormalSmooth)/44000))
            #print(bw)
            
            whno24 <- which(CGHTumor[,2]<=23)
            CGHTumor <- CGHTumor[whno24,]
            CGHTumor_large <- CGHTumor
            
            matched <- F
            if(nrow(CGHTumor) == nrow(CGHNormalSmooth)) if(sum(as.character(CGHTumor[,1]) != as.character(CGHNormalSmooth[,1]))==0) matched <- T
            
            if(!matched){
              mergeddata <- merge(CGHTumor_large,CGHNormalSmooth[,-c(2:n_end)],by=1,sort=F) #merges on probe name
              nnorm <- ncol(CGHNormalSmooth)
              ntum <- ncol(CGHTumor_large)
              CGHNormalSmooth <- mergeddata[,c(1:t_end,(ntum+1):(nnorm+ntum-n_end))]
              CGHTumor <- mergeddata[,1:ntum,drop=FALSE]
              
              #print("Different probe sets. Started probe matching")
              pll <- CGHTumor_large[,1]
              pls <- CGHTumor[,1]
              lps <- length(pls)
              lpl <- length(pll)
              ind <- 1:lps
              ind2 <- 1:lpl
              pls_ind <- data.frame(pls,ind)
              pll_ind <- data.frame(pll,ind2)
              #mergepr <- merge(pll_ind,pls_ind,by.x=1,by.y=1,all.x=TRUE,all.y=FALSE,sort=F) #merges on probe name
              mergepr <- match(pll_ind[,1],pls_ind[,1])
              #mergepord <- order(mergepr[,2])
              #indna <- mergepr[mergepord,3]
              indna <- mergepr
              ProbeBefore <- rep(NA,lpl)
              ProbeBefore[1] <- 1; pc <- 1
              for(i in 2:lpl){
                ind <- indna[i]
                if(!is.na(ind)) {pc <- ind}
                ProbeBefore[i] <- pc
              } 
              rindna <- rev(indna)
              ProbeAfter <- rep(NA,lpl)
              ProbeAfter[1] <- lps; pc <- lps
              for(i in 2:lpl){
                ind <- rindna[i]
                if(!is.na(ind)) {pc <- ind}
                ProbeAfter[i] <- pc
              }
              ProbeAfter <- rev(ProbeAfter)
              matchedProbes <- data.frame(pll,ProbeBefore,ProbeAfter)
              colnames(matchedProbes)<- c("probe","probeBefore","probeAfter")   
              probesNorm <- as.character(CGHNormalSmooth[,1])
              probesQuasi <- CGHTumor[,1]
              whichones <- match(probesNorm,probesQuasi)
              CGHNormalSmooth <- CGHNormalSmooth[!is.na(whichones),]
              whichoneind <- whichones[!is.na(whichones)]
            } #else {print("Same probe sets, matching not needed")}
            
            CGHann_large <- CGHTumor_large[,(1:t_end)]
            CGHTumor_large <- CGHTumor_large[,-(1:t_end),drop=FALSE]
            chromos_large <- CGHann_large[,2]
            
            
            print("Started pre-processing")
            chromos <- CGHNormalSmooth[,2]
            chrl <- as.character(sort(unique(chromos)))
            CGHann <- CGHTumor[,1:t_end]
            CGHNormalSmooth <- CGHNormalSmooth[,-(1:t_end)]
            CGHTumor <- CGHTumor[,-(1:t_end),drop=FALSE]
            
            nnorm <- ncol(CGHNormalSmooth)
            ntum <- ncol(CGHTumor)
            nboth <- nnorm+ntum
            
            CGHboth <- data.frame(CGHNormalSmooth,CGHTumor)
            
            #subtract and store median X-chromo  + global median normalization
            #whX <- which(chromos==23)
            #CGHX <- CGHboth[whX,]
            CGHnoX <- CGHboth
            ##CGHnoX <- CGHboth[-whX,]
            mediansnoX <- apply(CGHnoX,2,median)
            #mediansX <- apply(CGHX,2,median)
            #mediansXTumor <- mediansX[(nnorm+1):nboth]
            #CGHX <-  t(t(CGHX)-mediansX)
            #CGHboth[whX,] <- CGHX
            CGHboth <- t(t(CGHboth)-mediansnoX)
            CGHNormalSmooth <- CGHboth[,1:nnorm]
            CGHTumor <- CGHboth[,(nnorm+1):nboth,drop=FALSE]
            
            CGHTumor_large <- t(t(CGHTumor_large)-mediansnoX[(nnorm+1):nboth])
            
            #calculate moving mean
            print("Started computation of moving averages for tumor profiles")
            TumorsMovMean <- apply(CGHTumor,2,running,fun=mean,allow.fewer=TRUE,width=bw,align="left")  
            sdmean <- apply(TumorsMovMean,2,mad)
            
            Nopeaks <- function(prof, connorm=3){
              #prof <- NormalsSm0[,2];connorm<-3
              sdval <- mad(prof)
              wh0 <- which(abs(prof) >= connorm*sdval)
              prof[wh0] <- 0
              return(prof) 
            }
            
            NormalsSm <- apply(CGHNormalSmooth,2,Nopeaks,connorm = 3)
            
            print("Started robust estimation of sd for smoothed normal profiles")
            sdmeanNor <- apply(NormalsSm[seq(1,nrow(NormalsSm),by=bw),],2,running,fun=sd,allow.fewer=F,width=50,by=50) #skip bw probes, because these are smoothed values
            
            sdmeanNormed <- mean(apply(sdmeanNor,2,median))
            Vmeansig<- (sdmeanNormed)^2*nrow(NormalsSm)
            ridgepenal <- c(0,0.1,0.5,1,5,10,50,100,500,1000)*Vmeansig
            
            
            
            averMP <- function(i,mpi,allp){
              #mpi<-mp
              mp1 <- mpi[i,2]
              mp2 <- mpi[i,3]
              return((allp[mp1]+allp[mp2])/2)
            }
            
            enlarge <- function(allp=allpred,mp=matchedProbes){
              allpl <- sapply(1:nrow(mp),averMP,mpi=mp,allp=allp)
              return(allpl)
            }
            
            CorrectTumors_large_Fun <- function(tr,mat=matched,ridgereg=applyridgereg,ridgepen=ridgepenal,exclude=FALSE,sdm=sdmean,NormalsReg = NormalsSm,Tumors_large=CGHTumor_large,mp=matchedProbes,TumorsReg = CGHTumor,TumorsMM = TumorsMovMean,con=2.5,whichX=whX,medianX=mediansXTumor){
              #tr <-1;ridgereg=applyridgereg;ridgepen=ridgepenal;exclude=FALSE;sdm=sdmean;NormalsReg = NormalsSm;Tumors_large=CGHTumor_large;mp=matchedProbes;TumorsReg = CGHTumor;TumorsMM = TumorsMovMean;con=2.5;whichX=whX;medianX=mediansXTumor
              print(paste("Correcting tumor profile:",tr))
              if(exclude) {NormalsReg <- NormalsReg[,-tr]}
              whichdel <- which(abs(TumorsMM[,tr]) >= max(0.25,con*sdm[tr]))
              length(whichdel)
              if(length(whichdel)>0){NormalsRegDel <- NormalsReg[-whichdel,]} else {NormalsRegDel <- NormalsReg}
              TumorProf <- TumorsReg[,tr]
              if(length(whichdel)>0){TumorProfDel <- TumorProf[-whichdel]} else {TumorProfDel <- TumorProf}
              #NormalsReg <- Normals
              datatprof1 <- data.frame(y=TumorProfDel,NormalsRegDel)
              if(ridgereg) {
                reslqs <- lm.ridge(y ~ 0 + . , lambda=ridgepen,data = datatprof1)
                lambda_opt <- reslqs$kLW
              } else {
                lambda_opt<-0
              }
              reslqs <- lm.ridge(y ~ 0 + . , lambda=lambda_opt,data = datatprof1)
              scaled <-  reslqs$coef/reslqs$scales
              allpred <- as.matrix(NormalsReg)%*%scaled 
              if(matched){ 
                TP_large <- Tumors_large[,tr]
                TPcorrect_large <- TP_large - allpred
              } else {
                allpred_large <- enlarge(allpred,mp)
                TP_large <- Tumors_large[,tr]
                TPcorrect_large <- TP_large - allpred_large
              }   
              if(normalization != "median") {TPcorrect_large <- TPcorrect_large + mediansnoX[nnorm+tr]} #undo median normalization
              if(normalization == "mode") {
                denauto <- density(TPcorrect_large);
                mlauto <- denauto$x[which(denauto$y==max(denauto$y))] 
                TPcorrect_large <- TPcorrect_large - mlauto
              }
              return(TPcorrect_large)
            }
            
            
            corrected <- sapply(1:ntum,CorrectTumors_large_Fun,con=thr) 
            CGHcorrected_large <- data.frame(CGHann_large,corrected)
            colnames(CGHcorrected_large) <- c(colnames(CGHann_large),colnames(CGHTumor_large))        
            print("Finished")
            return(CGHcorrected_large)
          
      })

# EOF