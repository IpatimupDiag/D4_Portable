library(Biobase)
library(QDNAseq)

options(scipen = 999) #disabling scientific notation

#============================================================================
pjct <- snakemake@params[["projname"]]
RDS <- snakemake@input[["RDS"]]
profiles_dir <- snakemake@output[["table_out"]]
#============================================================================


makeSegmentationTable <- function(profile){
  #Get the profiles name:
  name <- profile$name

  #check bin size:
  bin_size =  (featureData(profile)$start[2] -1) / 1000

  #Status check:
  print(paste0("Sample: ",name))
  print(paste0("Bin size: ",bin_size,"Kbps"))

  #Get the all the used segment, and extrats their "coords" and "normalised log2 read count" values
  df_segmentationData <- data.frame(
    coords = rownames(assayData(profile)$segmented),
    segment_var = assayData(profile)$segmented)

  #Split Chr_Start_End
  split_coords <- strsplit(as.character(df_segmentationData$coords), split = ":")
  sub_df_segData <- do.call(rbind, split_coords) 
  #Split Start_End
  split_Start_End <- strsplit(as.character(sub_df_segData[,2]), split = "-")
  sub_sub_df_segData <- do.call(rbind, split_Start_End) 

  # Have to log2 segmentation data to be such as the profiles plot.
  #df_segmentationData <- data.frame(Chr = sub_df_segData[,1], Start_End = sub_df_segData[,2], Segment_var = log2(df_segmentationData[,2])) 
  df_segmentationData <- data.frame(
    Chr = sub_df_segData[,1], 
    Start = as.numeric(sub_sub_df_segData[,1]) , 
    End = as.numeric(sub_sub_df_segData[,2]) ,
    Segment_var = log2(df_segmentationData[,2]))

  #remove row with NA in Segment_var:
  df_segmentationData<-df_segmentationData[!is.na(df_segmentationData$Segment_var), ]

  # Merging intervals/Segments
  df_result <- data.frame(Chr = character(), Start = character(), End = character(), Segment_var = numeric(), stringsAsFactors = FALSE)

  i <- 1
  while (i <= nrow(df_segmentationData)) {
    chr <- df_segmentationData$Chr[i]
    seg_var <- df_segmentationData$Segment_var[i]
    start <- df_segmentationData$Start[i]
    end <- df_segmentationData$End[i]

    # Merge intervals while contiguous and have the same Chr & Segment_var
    while (i + 1 <= nrow(df_segmentationData) && df_segmentationData$Chr[i + 1] == chr && df_segmentationData$Segment_var[i + 1] == seg_var) {
      end <- df_segmentationData$End[i + 1]
      i <- i + 1
    }
  
    # Store the merged interval
    df_result <- rbind(df_result, data.frame(Chr = chr, Segment_Start = start, Segment_End = end, Segment_NormReadCount = seg_var))
  
    i <- i + 1
  }


  #Merge all bins with the same Segment_var in order, therefore defining the segments found in the profiles



  ###Save as a csv file
  ###file_name <- paste0(name,"-SegmentsTable-",bin_size)

  ###write.table(df_result,paste0(file_name,"kbp.tsv"),row.names=FALSE,sep="\t")
  #write.table(df_result,outputTable,row.names=FALSE,sep="\t")
  
  return(df_result)
}

data <- readRDS(RDS)

dir.create(profiles_dir)

for (i in 1:ncol(data)) {
  sn <- sampleNames(data)[i]

  tableOut_tsv<- paste(profiles_dir,paste(pjct, sn, "tsv", sep="."), sep="/") 
  print(tableOut_tsv)
  df_Seg_table <- makeSegmentationTable(data[,i])
  write.table(df_Seg_table,tableOut_tsv,row.names=FALSE,sep="\t")
}

#End:
print("All done!")
