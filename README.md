# ANext
Snakemake workflow to generate RNA-Seq assemblies and junctions

## How to run on High-Performance Computing (HPC) cluster
An example command is provided below for a SLURM scheduler
```bash
sbatch -c 1 --mem 5G -p ei-long -o out_ANext.%j.log -J ANext --wrap "source snakemake-5.4.0 && snakemake --latency-wait 120 --cluster-config /path/to/cluster.json --configfile /path/to/config.yaml --snakefile ANext.smk -p --jobs 100 --cluster \"sbatch -p {cluster.partition} -c {cluster.c} --mem {cluster.mem} -J {cluster.J} -o {cluster.o}\""
```
## Requirements
```
snakemake-5.4.0 - https://snakemake.readthedocs.io/en/stable/
hisat-2.1.0 - https://ccb.jhu.edu/software/hisat2/index.shtml
samtools-1.7 - http://www.htslib.org/download/
stringtie-1.3.3 - https://ccb.jhu.edu/software/stringtie/#install
scallop-0.10.2 - https://github.com/Kingsford-Group/scallop/releases
portcullis-1.1.2 - https://portcullis.readthedocs.io/en/latest/
```
## The config file runs the snakemake workflow, so edit it as needed
config file: config.yaml

## DAG
The below graph explains the basic workflow
![alt text](https://github.com/gemygk/ANext/blob/master/ANext.smk.svg)

