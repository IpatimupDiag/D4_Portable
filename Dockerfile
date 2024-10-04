##Set OS
FROM docker.io/continuumio/miniconda3:24.1.2-0

MAINTAINER Ze_Pedro

##Tries to solve the problem of snakemake not making dirs:
RUN mkdir /.Rcache; chmod a+rwX /.Rcache

#RUN mkdir folder

##Clone D4_Portable pipeline
#RUN git clone --single-branch --branch BETAv2.4 https://github.com/IpatimupDiag/D4_Portable/
RUN git clone https://github.com/IpatimupDiag/D4_Portable/

##Set path to work within the packages setup
WORKDIR /D4_Portable

##install packages from yaml-file
RUN conda config --set channel_priority flexible
RUN conda install -c conda-forge mamba
RUN mamba env update -n base --file env_D4_Portable.yml

##install non-conda R-dependencies
RUN Rscript r-dependencies.R

