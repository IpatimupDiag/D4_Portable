#!/usr/bin/env Rscript
##############################################################################################################
# script for QDNAseq profiles, gene representation over the CNA profile
# date: May 2025
##############################################################################################################

msg <- snakemake@params[["suppressMessages"]]
if (msg){
suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))
} else{
library(QDNAseq)
library(Biobase)
}

#---------------------------------------Functions:--------------------------------------------------
Find_Starting_Bin <- function(chromosome=22,profile){
	#Count to total of bins from the start util to the selected chr
	#And then convert the total o bins to the distance in the plot
	#Requires: the selected profile and the total of chromosome needed
	#Returns: a list with all starting postions of each chromossome 
	bin_size =  (featureData(profile)$start[2] -1) #/ 1000

	table_P <- featureData(profile)

	step <- 0
	mark_points_chr <- c()

	# start from 0
	for (i in 0:(chromosome -1)){
		if (i != 0){
			#Just the bins from the selected chr:
			sub_table_chr <- table_P[table_P$chromosome==i,]

			#Getting the end coords of the last bin:
			endPos_last_bin <- sub_table_chr[nrow(sub_table_chr)]$end

			#Sums the last marked step with the end coords of the last bin
			#Saves the last step from the he end coords of the last bin
			step <- step + endPos_last_bin
		}
		#Then had this new marker where the chr ended
		mark_points_chr <- append(mark_points_chr, step)
	}
	return(mark_points_chr)
}

ByEach <- function(profile,df_OncoGenes){
	# By each autosome (1 - 22) will get the staring postion in the QDNAseq plot

	#Getting a list of the start positions where the first bin of each chr starts in the plot:
	list_start_binPos <- Find_Starting_Bin(22,profile)

	for (chr in seq(22)){
		# Now get the Onco genes that belong to the same chromosome:
		df_matching_chr_genes <- df_OncoGenes[df_OncoGenes$Chromosome == chr,]

		#Getting the position where the first bin of the selected chr stars in the plot:
		start_binPos <- list_start_binPos[chr]
		
		#
		# Check is there is a onco gene in the selected chr:
		if (dim(df_matching_chr_genes)[1] != 0){

			# One by one Onco gene, from the same chr, is marked in the plot:
			for (n_gene in 1:nrow(df_matching_chr_genes) ){
				gene_name <- df_matching_chr_genes[n_gene,]$Genes
				gene_str <- df_matching_chr_genes[n_gene,]$Start
				gene_end <- df_matching_chr_genes[n_gene,]$End
				gene_col <- df_matching_chr_genes[n_gene,]$info
				gene_level <- df_matching_chr_genes[n_gene,]$level
			
				#Now mark the gene in the plot:
				GeneBookmark(gene_name,gene_str,gene_end,gene_col,gene_level,start_binPos)
			}
		}
	}
}



GeneBookmark <- function(gene_name,gene_str,gene_end,gene_col,gene_level,start_binPos){
	#' Will mark the selected OncoGene to the existing QDNAseq plot as a strait vertical line

	# marking the start and end of the gene in the plot
	gene_start_plot <- start_binPos+ gene_str # marks the start of the chromosome in the plot
	gene_end_plot <- start_binPos+ gene_end # marks the end of the chromosome in the plot

	# draw gene in plot
	ymin <- 1.5 #*gene_level # (multipy by gene_level, by 2 or 1)
	ymax <- 2 + (gene_level)
	rect(gene_start_plot, ymin, gene_end_plot, ymax, col = gene_col, border = gene_col)

	# add label for gene name
	gene_name <- gene_name
	gene_label_x <- (gene_start_plot + gene_end_plot) / 2
	gene_label_y <- ymax + 0.1 * (ymax - ymin)
	text(gene_label_x, gene_label_y, labels = gene_name, col = gene_col, cex = 0.85, pos = 3,srt=90)

}

#---------------------------------------\Functions:--------------------------------------------------

#----------------------------------INPUT/OUTPUT-------------------------------------------------------------------------
# consider if is runinng in a snakemake pipeline.
if (exists("snakemake")){
	input_QDNAseq <- snakemake@input[["RDS"]]
	profiles_dir <- snakemake@output[["profiles_path"]]
	file_OncoGenesList <- snakemake@params[["genes_list"]]
}else{
	args = commandArgs(trailingOnly=TRUE)
	input_QDNAseq <- args[1] # INPUT # ex: "./test_Case/100kbp-segmented.rds" # Getting the copy number profiles:
	profiles_dir <- args[2] # OUTPUT # ex: "./.../GeneRuler_Info/"
	file_OncoGenesList <- "data_GeneRuler/GenesRuler_hg19-Ensembl_allGenes.csv"
}
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Read the OncoGenes list from a csv file:
df_OncoGenes <- read.csv(file_OncoGenesList, header=1,sep=",")

# Read the copy number profiles:
QRC <- readRDS(input_QDNAseq)
#------------------------------------------------------------------------------------
# Checking the bin size:
if (grepl("100k",input_QDNAseq, fixed=TRUE)){
	bin_num <- 100000 # 100K
	bin_size <- "100k"
} else if (grepl("1000k",input_QDNAseq, fixed=TRUE)) {
	bin_num <- 1000000 #1000k
	#bin_num <- 100000
	bin_size <- "1000K"
	#print("----- GOINING 1000K BINS !!!!!")
} else {
	print("Wrong file !!! No bin size in the file's name")
	print("The input file should be something like 100kbp-segmented.rds or 1000kbp-segmented.rds")
}

# Make the missing directory
dir.create(profiles_dir)
#------------------------------------------------------------------------------------------------------------
# Working by each copy number profile avaliable at the QDNAseq file:
total_profiles <- length(sampleNames(QRC))

for (i in seq(total_profiles)) {
	current_profile <-QRC[,i]

	# Prepare to save plot:
	profile_name <- sampleNames(current_profile)

	print(paste0("Ploting ",profile_name))

	#png(paste(profiles_dir,paste(bin_size,"bps_",profile_name, ".OncoGeneRuler.TESTMar25-TEST100.png", sep=""), sep=""), width=1690, height=480)
	png(paste(profiles_dir,paste("/",bin_size,"bps_Bin_",profile_name, ".OncoGeneRuler.Mar25-allGenes.png", sep="")
		, sep="")
	, width=600, height=300, units='mm', res=250)

	# Make the plot:
	plot(current_profile,  ylim=c(-2,5))

	# Check for matching Onco genes and marks then to the profile:
	ByEach(current_profile,df_OncoGenes)

	#
	dev.off()
}


