### R code from vignette source 'Clonality.Rnw'

###################################################
### code chunk number 1: read data
###################################################
 library(Clonality)
set.seed(100)
chrom<-paste("chr",rep(c(1:22),each=100),"p",sep="")
chrom[nchar(chrom)==5]<-paste("chr0",substr(chrom[nchar(chrom)==5] ,4,5),sep="")
maploc<- rep(c(1:100),22)
data<-NULL
for (pt in 1:9)  #first 9 patients have independent tumors
{
tumor1<-tumor2<- NULL
mean1<- rnorm(22) 
mean2<- rnorm(22)
for (chr in 1:22)
{ 
  r<-runif(2) 
if (r[1]<=0.5) tumor1<-c(tumor1,rep(0,100))   
  else if   (r[1]>0.7)  tumor1<-c(tumor1,rep(mean1[chr],100))
  else  { i<-sort(sample(1:100,2))
        tumor1<-c(tumor1,mean1[chr]*c(rep(0,  i[1]),rep(1, i[2]-i[1]), rep(0,  100-i[2])))
        }
if (r[2]<=0.5) tumor2<-c(tumor2,rep(0,100))
  else if   (r[2]>0.7)  tumor2<-c(tumor2,rep(mean2[chr],100))
  else   {i<-sort(sample(1:100,2))
       tumor2<-c(tumor2,mean2[chr]*c(rep(0,  i[1]),rep(1, i[2]-i[1]), rep(0,  100-i[2])))
         }
}
data<-cbind(data,tumor1,tumor2)
}

#last patient has identical profiles
tumor1<- NULL
mean1<- rnorm(22) 
for (chr in 1:22)
{ 
  r<-runif(1) 
if (r<=0.4) tumor1<-c(tumor1,rep(0,100))   
  else if   (r>0.6)  tumor1<-c(tumor1,rep(mean1[chr],100))
  else  { i<-sort(sample(1:100,2))
        tumor1<-c(tumor1,mean1[chr]*c(rep(0,  i[1]),rep(1, i[2]-i[1]), rep(0,  100-i[2])))
        }

}
data<-cbind(data,tumor1,tumor1)

data<-data+matrix(rnorm( 44000,mean=0,sd=0.4) ,nrow=2200,ncol=20)


samnms<-paste("pt",rep(1:10,each=2),rep(1:2,10),sep=".")



###################################################
### code chunk number 2: Clonality.Rnw:105-106
###################################################
dim(data)


###################################################
### code chunk number 3: CNA
###################################################

dataCNA<-CNA(data,chrom=chrom,maploc=maploc,sampleid=samnms)
as.matrix(dataCNA)[1:5,1:10]



###################################################
### code chunk number 4: averaging
###################################################
dataAve<- ave.adj.probes(dataCNA,2)


###################################################
### code chunk number 5: Clonality.Rnw:131-132
###################################################
ptlist<- paste("pt",rep(1:10,each=2),sep=".")


###################################################
### code chunk number 6: main clonality analysis
###################################################
results<-clonality.analysis(dataCNA, ptlist,  pfreq = NULL, refdata = NULL, nmad = 1,  reference = TRUE, allpairs = TRUE)


###################################################
### code chunk number 7: Clonality.Rnw:141-142
###################################################
results$LR


###################################################
### code chunk number 8: genomewidePlots
###################################################
genomewidePlots(results$OneStepSeg, results$ChromClass,ptlist , c("pt.10.1", "pt.10.2"),results$LR,  plot.as.in.analysis = TRUE) 


###################################################
### code chunk number 9: Clonality.Rnw:155-156
###################################################
chromosomePlots(results$OneStepSeg, ptlist,ptname="pt.10",nmad=1)


###################################################
### code chunk number 10: Clonality.Rnw:161-162
###################################################
histogramPlot(results$LR[,4], results$refLR[,4])


###################################################
### code chunk number 11: LOH data generation
###################################################
set.seed(25)
LOHtable<-cbind(1:20,matrix(sample(c(0,1,2),20*20,replace=TRUE),20))
LOHtable[,3]<-LOHtable[,2]
LOHtable[1,3]<-0


###################################################
### code chunk number 12: LOH data analysis
###################################################
LOHtable[,1:5]
LOHclonality(LOHtable,rep(1:10,each=2),pfreq=NULL,noloh=0,loh1=1,loh2=2)


###################################################
### code chunk number 13: sessionInfo
###################################################
sessionInfo()


