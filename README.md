# ANext
Snakemake workflow to generate RNA-Seq assemblies and junctions

## DAG
The below graph explains the basic workflow
![alt text](https://github.com/gemygk/ANext/blob/master/ANext.smk.svg)

## Requirements
```
snakemake-5.4.0 - https://snakemake.readthedocs.io/en/stable/
hisat-2.1.0 - https://ccb.jhu.edu/software/hisat2/index.shtml
samtools-1.7 - http://www.htslib.org/download/
stringtie-1.3.3 - https://ccb.jhu.edu/software/stringtie/#install
scallop-0.10.2 - https://github.com/Kingsford-Group/scallop/releases
portcullis-1.1.2 - https://portcullis.readthedocs.io/en/latest/
```
## Clone workflow into working directory
```bash
git clone https://github.com/gemygk/ANext.git path/to/workdir
cd path/to/workdir
```

## Test the workflow
Once you have working snakemake installed, please test the workflow
```bash
cd path/to/workdir
snakemake --snakefile ANext.smk --configfile config.yaml -np
```
The above should print the following:
```
Building DAG of jobs...
...
...
...
Job counts:
	count	jobs
	1	all
	3	hisat2
	1	hisat2_build
	1	portcullis
	3	scallop
	3	stringtie
	12
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

## The config file runs the workflow, so edit it as needed
config file: [config.yaml](https://github.com/gemygk/ANext/blob/master/config.yaml)

## Run on High-Performance Computing (HPC) cluster
An example command is provided below for [SLURM](https://slurm.schedmd.com/) scheduler
```bash
cd path/to/workdir
sbatch -c 1 --mem 5G -p ei-long -o out_ANext.%j.log -J ANext --wrap "source snakemake-5.4.0 && snakemake --latency-wait 120 --cluster-config cluster.json --configfile config.yaml --snakefile ANext.smk -p --jobs 100 --cluster \"sbatch -p {cluster.partition} -c {cluster.c} --mem {cluster.mem} -J {cluster.J} -o {cluster.o}\""
```
Please modify the cluster configuration file [cluster.json](https://github.com/gemygk/ANext/blob/master/cluster.json) as needed to work on any other HPC scheduler (like PBS Pro, LSF etc)


