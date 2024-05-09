#Set OS
FROM docker.io/continuumio/miniconda3:24.1.2-0

MAINTAINER Ze_Pedro

##RUN mkdir folder
RUN git clone https://github.com/jppmatos/D4_Portable/

#Set path to work within the packages setup
WORKDIR /D4_Portable

#install packages from yaml-file
RUN conda config --set channel_priority flexible
RUN conda install -c conda-forge mamba
RUN mamba env update -n base --file env_D4_BETA.yml

#install non-conda R-dependencies
RUN Rscript r-dependencies.R

