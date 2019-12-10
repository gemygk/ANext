#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Snakemake workflow to generate RNA-Seq assemblies and junctions

21/02/2019, 10:58:11
"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__copyright__ = "Copyright 2019"
__license__ = "GNU General Public License v3.0"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "Gemy.Kaithakottil@gmail.com, gemygk@gmail.com"
__status__ = "Production"
__version__ = "0.01"

# import modules
import os
import sys
from glob import glob
import subprocess


# Request min version of snakemake
from snakemake.utils import min_version
min_version("5.4.0")

# declare variables
cwd = os.getcwd()

# get genome
genome = os.path.abspath(config["genome"])
if not os.path.exists(genome):
    print ("ERROR: The genome file cannot be accessed - '{0}'".format(genome))
    sys.exit()
genome_base = os.path.basename(genome)
# check if genome is gzipped and exit if true
gzipped = None
gzipped = True if genome.lower().endswith(".gz") else False
if gzipped:
    print("ERROR: '{0}' file should not be gzipped, please unzip the file and rerun".format(genome))
    sys.exit()

# check if genome is over 2^31-1, if yes, create CSI index
cmd = "awk 'BEGIN {total=0} {if($0 !~ /^>/) {total+=length}} END{print total}' " +  genome
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
(result,error) = p.communicate()
exit_code = p.returncode
if exit_code:
    raise subprocess.CalledProcessError(exit_code,cmd)
result = int(result) # convert str to int
print("The worked out genome length is:{}".format(result))
csi_index = None
if result < (2**31)-1:
    # print ("INFO:Normal bai index {} < 2**31-1".format(result))
    csi_index = False
else:
    # print ("INFO:Create csi index {} > 2**31-1".format(result))
    csi_index = True

# get tools
HISAT_VERSION = None
STRINGTIE_VERSION = None
SCALLOP_VERSION = None
PORTCULLIS_VERSION = None
load = config["load"]
for tool, source in load.items():
    if tool in ("hisat2"):
        source_dump, HISAT_VERSION = source.split(" ")
    if tool in ("stringtie"):
        source_dump, STRINGTIE_VERSION = source.split(" ")
    if tool in ("scallop"):
        source_dump, SCALLOP_VERSION = source.split(" ")
    if tool in ("portcullis"):
        source_dump, PORTCULLIS_VERSION = source.split(" ")

# get output folder
OUTPUT = os.path.abspath(config["output"])

# get sample information
# work out if samples are gzipped or not
rna_seq_samples = []
read_dir = os.path.join(OUTPUT,"data","reads")
all_samples = config["samples"]
gzipped = None
for sample_name, reads in all_samples.items():
    rna_seq_samples.append(sample_name)
    reads_dict = {}
    reads_dict["R1"] = reads[0]
    reads_dict["R2"] = reads[1]
    for pair, path in reads_dict.items():
        r_path = os.path.abspath(path)
        gzipped = True if r_path.lower().endswith(".gz") else False
        if not os.path.exists(r_path):
            print ("ERROR: The pair '{0}' of read '{1}' for sample '{2}' cannot be accessed".format(pair,r_path,sample_name))
            sys.exit()
        if not os.path.exists(read_dir):
            os.makedirs(read_dir)
        new_read_name_list = (sample_name,pair,"fastq.gz") if gzipped else (sample_name,pair,"fastq")
        new_read_name = ".".join(new_read_name_list)
        cmd = "cd " + read_dir + " && ln -s " + r_path + " " + new_read_name
        if not os.path.exists(os.path.join(read_dir,new_read_name)):
            os.system(cmd)

# get strandedness
strandedness = config["strandedness"]
if strandedness not in ("fr-firststrand","fr-secondstrand", "unstranded"):
    print ("ERROR: The strand option provided '{0}' is not in vaild format. The valid formats are fr-firststrand, fr-secondstrand, OR unstranded".format(strandedness))
    sys.exit()

hisat2_strand = None # tool help here : https://ccb.jhu.edu/software/hisat2/manual.shtml
stringtie_strand = None # tool help here : http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
scallop_strand = None # tool help here : https://github.com/Kingsford-Group/scallop
portcullis_strand = None # tool help here : https://portcullis.readthedocs.io/en/latest/using.html
if strandedness in ("fr-firststrand"):
    hisat2_strand = "RF"
    stringtie_strand = "--rf"
    scallop_strand = "first"
    portcullis_strand = "firststrand"
elif strandedness in ("fr-secondstrand"):
    hisat2_strand = "FR"
    stringtie_strand = "--fr"
    scallop_strand = "second"
    portcullis_strand = "secondstrand"
else:
    hisat2_strand = "unstranded"
    stringtie_strand = "unstranded"
    scallop_strand = "unstranded"
    portcullis_strand = "unstranded"

# create logs folder
# need to find a proper fix for this as mentioned in the issue below, but for now using a quick fix
# # https://bitbucket.org/snakemake/snakemake/issues/838/how-to-create-output-folders-for-slurm-log#comment-45348663
# for rule in workflow._rules:
#     getattr(rules, rule)._set_params_item(config["my_slurm_dir"], name = "slurm_dir")
# os.makedirs(config["my_slurm_dir"], exist_ok = True)
cluster_logs_dir = os.path.join(cwd,"logs","cluster")
if not os.path.exists(cluster_logs_dir):
    os.makedirs(cluster_logs_dir)


#######################
# RULES STARTS HERE
#######################

shell.prefix("set -eo pipefail; ")

rule all:
    input:
        expand(os.path.join(OUTPUT,HISAT_VERSION,"index",genome_base + ".hisat2-build.done")),
        expand(os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.{ext}"),sample=rna_seq_samples,ext=['bam', 'bam.csi']) if csi_index else expand(os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.{ext}"),sample=rna_seq_samples,ext=['bam', 'bam.bai']),
        expand(os.path.join(OUTPUT,STRINGTIE_VERSION,"{sample}","stringtie_{sample}.transcripts.gtf"),sample=rna_seq_samples),
        expand(os.path.join(OUTPUT,SCALLOP_VERSION,"{sample}","scallop_{sample}.transcripts.gtf"),sample=rna_seq_samples),
        expand(os.path.join(OUTPUT,PORTCULLIS_VERSION,"portcullis_out","3-filt","portcullis_filtered.{status}.junctions.{ext}"),
            status = ['pass', 'fail'], ext = ['tab', 'bed'])


#######################
# ALIGNMENT & ASSEMBLY WORKFLOW
#######################

# run hisat2-build
# ------------------
rule hisat2_build:
    input:
        genome = genome
    output:
        expand(os.path.join(OUTPUT,HISAT_VERSION,"index",genome_base + ".hisat2-build.done")),
    log:
        cwd = expand(os.path.join(OUTPUT,HISAT_VERSION,"index","hisat2_build.log"))
    threads: 4
    params:
        source = config["load"]["hisat2"],
        cwd = expand(os.path.join(OUTPUT,HISAT_VERSION,"index"))
    shell:
        "(set +u"
        + " && cd {params.cwd}"
        + " && ln -s {input.genome}"
        + " && {params.source} && /usr/bin/time -v hisat2-build -p {threads} {genome_base} {genome_base}" \
        + " && touch {output}) 2> {log}"

# run hisat2
# -----------
# get input pattern - really helpful biostars link
# https://www.biostars.org/p/296020/#296132
def get_r1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    return glob(os.path.join(read_dir,wildcards.sample + '.R1.fastq*'))

def get_r2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample* e.g.
    return glob(os.path.join(read_dir,wildcards.sample + '.R2.fastq*'))

rule hisat2:
    input:
        lambda wildcards: config["samples"][wildcards.sample],
        r1= get_r1,
        r2= get_r2,
        idx = os.path.join(OUTPUT,HISAT_VERSION,"index",genome_base + ".hisat2-build.done")
    output:
        unsorted_bam = os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.dta.unsorted.bam"),
        bam = os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.bam"),
        bai = os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.bam.csi") if csi_index else os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.bam.bai")
    log:
        os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.hisat2.log")
    params:
        cwd = os.path.join(OUTPUT,HISAT_VERSION,"{sample}"),
        idx = expand(os.path.join(OUTPUT,HISAT_VERSION,"index","{genome}"),genome=genome_base),
        source_hisat2 = config["load"]["hisat2"],
        source_samtools = config["load"]["samtools"],
        extra = config["load_parameters"]["hisat2"] + " --rna-strandness " + hisat2_strand if hisat2_strand not in ("unstranded") else config["load_parameters"]["hisat2"],
        csi_index_status = "-c" if csi_index else ""
    threads: 4
    shell:
        "(set +u" \
        + " && cd {params.cwd} "
        + " && {params.source_hisat2} "
        + " && {params.source_samtools} "
        + " && /usr/bin/time -v hisat2 {params.extra} -x {params.idx} -p {threads} -1 {input.r1} -2 {input.r2} "
        + " | /usr/bin/time -v samtools view -b -@ {threads} - > {output.unsorted_bam} "
        + " && /usr/bin/time -v samtools sort --threads {threads} {output.unsorted_bam} > {output.bam} "
        + " && /usr/bin/time -v samtools index -@ {threads} {params.csi_index_status} {output.bam}) "
        + " 2> {log}"


# run stringtie
# -----------
rule stringtie:
    input:
        bam = os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.bam")
    output:
        os.path.join(OUTPUT,STRINGTIE_VERSION,"{sample}","stringtie_{sample}.transcripts.gtf")
    log:
        os.path.join(OUTPUT,STRINGTIE_VERSION,"{sample}","stringtie_{sample}.log")
    params:
        cwd = os.path.join(OUTPUT,STRINGTIE_VERSION,"{sample}"),
        source = config["load"]["stringtie"],
        extra = config["load_parameters"]["stringtie"] + " " + stringtie_strand if stringtie_strand not in ("unstranded") else config["load_parameters"]["stringtie"]
    threads: 4
    shell:
        "(set +u"
        + " && cd {params.cwd} "
        + " && {params.source} "
        + " && /usr/bin/time -v stringtie {input.bam} {params.extra} -p {threads} -o {output}) 2> {log}"


# run scallop
# -----------
rule scallop:
    input:
        bam = os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.bam")
    output:
        os.path.join(OUTPUT,SCALLOP_VERSION,"{sample}","scallop_{sample}.transcripts.gtf")
    log:
        os.path.join(OUTPUT,SCALLOP_VERSION,"{sample}","scallop_{sample}.log")
    params:
        cwd = os.path.join(OUTPUT,SCALLOP_VERSION,"{sample}"),
        source = config["load"]["scallop"],
        extra = config["load_parameters"]["scallop"] + " --library_type " + scallop_strand if scallop_strand not in ("unstranded") else config["load_parameters"]["scallop"]
    threads: 4
    shell:
        "(set +u"
        + " && cd {params.cwd} "
        + " && {params.source} "
        + " && /usr/bin/time -v scallop -i {input.bam} -o {output} {params.extra}) 2> {log}"


# run portcullis
# -----------
rule portcullis:
    input:
        genome = genome,
        bam = expand(os.path.join(OUTPUT,HISAT_VERSION,"{sample}","{sample}.bam"),sample=rna_seq_samples),
    output:
        expand(os.path.join(OUTPUT,PORTCULLIS_VERSION,"portcullis_out","3-filt","portcullis_filtered.{status}.junctions.{ext}"),
        status = ['pass', 'fail'], ext = ['tab', 'bed']),
    log:
        os.path.join(OUTPUT,PORTCULLIS_VERSION,"portcullis_out","portcullis.log")
    params:
        cwd = os.path.join(OUTPUT,PORTCULLIS_VERSION),
        source = config["load"]["portcullis"],
        extra = config["load_parameters"]["portcullis"] + " --output portcullis_out --strandedness " + portcullis_strand,
        csi_index_status = "--use_csi" if csi_index else ""
    threads: 8
    shell:
        "(set +u"
        + " && cd {params.cwd} "
        + " && {params.source} "
        + " && /usr/bin/time -v portcullis full --threads {threads} {params.csi_index_status} {params.extra} {input.genome} {input.bam}) 2> {log}"

#######################
# END
#######################
