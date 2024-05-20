# D4_portable
D4 pipeline, for Copy Number Profiles analysis

#### :warning: This pipeline is stil in developement, changes may occur in the pipeline and in the repository itself! :construction:


### Run within Podman:
```bash
podman run -v /shared/input_dir:/mnt/shared/:z docker.io/ipatimupdiag/d4portable snakemake --cores 2
```


##### Run in interactive mode:
```bash
podman run -it -v /shared/input_dir:/mnt/shared/:z docker.io/ipatimupdiag/d4portable bash
```
##### Run specific branch:
```bash
podman run -v /shared/input_dir:/mnt/shared/:z docker.io/ipatimupdiag/d4portable:BETAv2.5 snakemake --cores 2 --configfile /mnt/shared/config.yaml
```

##### Run with custom config file:
```bash
podman run -v /shared/input_dir:/mnt/shared/:z docker.io/ipatimupdiag/d4portable snakemake --cores 2 --configfile /mnt/shared/config.yaml 
```
