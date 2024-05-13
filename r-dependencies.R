#!/usr/bin/env Rscript
##############################################################################

#Author: Erik Bosch
#date: Jan 20221

#edited 
#Author: Jose Pedro Matos
#date: Apr 2024

#This script is for installing remote packages to be used in the R-environment for snakemake
##############################################################################

#install.packages("devtools", version="2.4.5",repos="http://cran.us.r-project.org")
#options(unzip = 'internal')

#if (!require("BiocManager", quietly = TRUE))
#            install.packages("BiocManager",repos="http://cran.us.r-project.org")
#BiocManager::install("DNAcopy", version="1.72.3")

install.packages("renv", repos='http://cran.us.r-project.org')

options(renv.settings.bioconductor.version = "3.16")
renv::install("bioc::DNAcopy@1.72.3")

# install WECCA from Github
devtools::install_github("tgac-vumc/WECCA",force=TRUE)

# install QDNAseq.dev
devtools::install_github("tgac-vumc/QDNAseq.dev", ref="clonality",force=TRUE)
#devtools::install_github("tgac-vumc/QDNAseq.dev@clonality", force=TRUE)
#devtools::install_github("tgac-vumc/QDNAseq.dev", ref="clonality",force=TRUE)

# install GGHtest from source
fn = 'http://www.few.vu.nl/~mavdwiel/CGHtest/CGHtest_1.1.tar.gz'
if (!file.exists('CGHtest_1.1.tar.gz')){download.file(fn, destfile="CGHtest_1.1.tar.gz")} else {untar("CGHtest_1.1.tar.gz",list=TRUE)} # download or show content ofalready downloaded

install.packages("CGHtest_1.1.tar.gz", source=NULL, depenencies=TRUE)
