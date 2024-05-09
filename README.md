# D4_portable
D4 pipeline, for Copy Number Profiles analysis

### This pipeline is stil in developement, changes may occur in the pipeline and in the repository itself!


##### Podman:

podman run -it -v /shared/input_dir:/mnt/shared/:Z docker.io/jpmatos/test:0.5.D4_Portable bash

podman run -u root -v /shared/input_dir:/mnt/shared/:Z docker.io/jpmatos/test:0.5.D4_Portable snakemake --cores 1 --configfile /mnt/shared/config.yaml 
