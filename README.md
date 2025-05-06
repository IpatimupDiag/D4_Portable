# D4_portable
A clonality identification tool for lung cancer lesions, that runs thru Copy Number alteration Profiles analysis.

## About
Lung cancer is the leading cause of cancer-related deaths worldwide, often presenting with multiple pulmonary lesions. Synchronous or metachronous lesions complicate treatment decisions, especially when standard diagnostics cannot define their relationships. Distinguishing between metastatic (clonal) and new primary (non-clonal) tumors is crucial for patient management. 

This pipeline implements molecular diagnostics method using DNA copy number alteration (CNA) profiles form shallow whole-genome (sWGS), Illumina and Ion Torrent technology, to determine the clonality of multiple lung cancer lesions, assisting clinical decision-making.



This project was performed under collaboration between the Institute of Molecular Pathology and Immunology of the University of Porto (IPATIMUP) and Amsterdam UMC's Tumor Genome Analysis Core (TGAC).
The current was based of TGAC's clonality pipeline and QDNAseq.snakemake pipeline.


### Requirements:
Just have podman (version **3.4.4** recommended) installed if you what to run it in container mode. Otherwise, have conda and mamba installed.



### How to run it!

#### Run within Podman:
```bash
podman run -v /<home_dir>/shared:/mnt/shared/:z docker.io/ipatimupdiag/d4portable snakemake --cores 2
```
#### or
```bash
podman run -v /<home_dir>/projData:/mnt/shared/input_dir/:z -v /<home_dir>/projReport:/mnt/shared/output_dir/:z  docker.io/ipatimupdiag/d4portable snakemake --cores 2
```


**input:**
```
/<home_dir>/shared/input_dir/bam/
```

**output:**
```
/<home_dir>/shared/output_dir/resultReport/
```


##### Run in interactive mode:
```bash
podman run -it -v /<home_dir>/shared:/mnt/shared/:z docker.io/ipatimupdiag/d4portable bash
```
##### Run specific branch:
```bash
podman run -v /<home_dir>/shared:/mnt/shared/:z docker.io/ipatimupdiag/d4portable:BETAv2.5 snakemake --cores 2 --configfile /mnt/shared/config.yaml
```

##### Run with custom config file:
```bash
podman run -v /<home_dir>/shared:/mnt/shared/:z docker.io/ipatimupdiag/d4portable snakemake --cores 2 --configfile /mnt/shared/config.yaml 
```

#### Sugestion of a bash script to run D4 Portable
```bash
#!/bin/sh

inputDir=$1':/mnt/shared/input_dir/bam/:z'
outputDir=$2':/mnt/shared/output_dir/:z'

podman run 
        --rm \
        -v $inputDir \
        -v $outputDir \
        docker.io/ipatimupdiag/d4portable \
        snakemake --cores 2
```
```bash
$ sh run_D4_Portable.sh /{your_path}/{input_dir} /{your_path}/{output_dir}/
```

#### Run within conda:
```bash
conda config --set channel_priority flexible
```
```bash
conda install -c conda-forge mamba
```
```bash
mamba env update -n D4_Portable --file env_D4_Portable.yml
```
```bash
conda activate D4_Portable
```
```bash
Rscript r-dependencies.R
```
###### You will have to update **config.yaml**, *path*: *dir_bam* and *dir_out* !
```bash
snakemake --cores 2 --configfile config.yaml
```


