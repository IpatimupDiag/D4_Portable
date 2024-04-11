#Set OS
FROM docker.io/continuumio/miniconda3:24.1.2-0

MAINTAINER Ze_Pedro

#RUN apt-get update -y

#RUN apt-get --qq -y install curl tar bzip2 git zip wget


#Pull the Pipeline
WORKDIR /

#RUN mkdir folder

RUN git clone https://github.com/tgac-vumc/QDNAseq.snakemake/ 

RUN ls

RUN cd QDNAseq.snakemake

#install packages from yaml-file
RUN conda install --file environment.yaml -y

#install non-conda R-dependencies
RUN Rscript r-dependencies.R


RUN conda list

RUN echo ''

RUN echo 'hello hello'

