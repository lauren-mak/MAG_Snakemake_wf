
from os.path import basename, exists
import pandas as pd


# From https://nf-co.re/mag/2.1.0/output 
rule porechop:
    input:
        join(FASTQ_DIR, "{run}.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "porechop/{run}.fastq.gz"),
    threads: workflow.cores,
    shell:
        """
        porechop -i {input} -t {threads} -o {output}
        """


# Parameters mostly from https://nf-co.re/mag/2.1.0/output 
rule filtlong:
    input:
        join(DATA_DIR, preprocessing_dir, "porechop/{run}.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
    params:
        short_qc="1000",
        long_len_wt="10",
    threads: workflow.cores,
    shell:
        """
        filtlong --min_length {params.short_qc} --keep_percent 90 --length_weight {params.long_len_wt} --trim --target_bases 500000000 {input} | gzip > {output}
        """


# rule nanoplot:
#     input:
#         join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
#     output:
#         TODO,
#     params:
#         outdir=join(DATA_DIR, preprocessing_dir, "postprocessing_qc"),
#     shell:
#         """
#         NanoPlot -o {params.outdir} --loglength -f pdf --N50 --legacy kde --fastq {input}
#         """


rule count_reads:
    input:
        join(DATA_DIR, preprocessing_dir, "singlerun/{run}_1.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "singlerun/{run}_readcount.csv"),        
    params:
        sample_name="{run}",
    shell:
        """
        num_reads=`cat {input} | awk 'NR%4==2{c++; l+=length($0)} END {print c","l}'`
        echo {sample_name},${num_reads} > {output}
        """
        # Adapted from https://wiki.itap.purdue.edu/display/CGSB/fastq%3A+count+reads+in+fastq+or+gzipped+fastq+file


def est_mem_flye(wildcards): # Unnecessary, for now
    info = open(join('data/00_preprocessing/singlerun', wildcards.run + '_readcount.csv'),"r").readline().split(',')
    num_reads = int(info[2]) # sample_name,num_reads,num_bases
    base_mem = (int)((0.00000183542567428533 * num_reads - 8.01103718491264) * 1.1) # [-] Different?
    attempt_f = join(DATA_DIR, assembly_dir,'singlerun', wildcards.run + '_attempt.txt') # Need this because wildcards.attempt is inaccessible?
    attempt_c = 0  
    if exists(attempt_f):
        attempt_c = (int)(open(attempt_f, 'r').readlines()[0].strip()) + 1
    print(attempt_c, file = open(attempt_f, 'w'))
    return max(max(base_mem, 0) + 100 * attempt_c, 50) # [-] Different?

rule metaflye:
    input:
        join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
        # join(DATA_DIR, preprocessing_dir, "singlerun/{run}_readcount.txt"),        
    output:
        join(DATA_DIR, assembly_dir, "flye_unpolished/singlerun/{run}/assembly.fasta"),
    threads: workflow.cores,
    params:
        outdir=join(DATA_DIR, assembly_dir, "flye_unpolished/singlerun/{run}"),
    shell:
        """
        if [ -f "{params.outdir}/flye.log" ]; then
            flye --resume -o {params.outdir} # Resume from the last previously completed step. Don't have to specify
        else
            flye --nano-raw {input} -o {params.outdir} -t {threads} --meta
        fi
        """


rule polish_mapping:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
        asm=join(DATA_DIR, assembly_dir, "flye_unpolished/singlerun/{run}/assembly.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{run}.bam"),
    threads: workflow.cores,
    params:
        sample_name="{run}",
        outdir=join(DATA_DIR, assembly_dir, "prelim_polished/singlerun"),
    shell:
        """
        minimap2 -ax map-ont {input.asm} {input.fq} | samtools view -bS - > {params.outdir}/{params.sample_name}_tmp.bam
        samtools sort -@ {threads} -o {output}
        rm {params.outdir}/{params.sample_name}_tmp.bam
        """

rule racon:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
        mmp=join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{run}.bam"),
        asm=join(DATA_DIR, assembly_dir, "flye_unpolished/singlerun/{run}/assembly.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{run}.racon.fasta"),
    threads: workflow.cores,
    shell:
        """
        racon --t {threads} {input.fq} {input.mmp} {input.asm} > {output}
        """


rule medaka:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
        cn=join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{run}.racon.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{run}.medaka.fasta"),
    threads: workflow.cores,
    params:
        model="r941_prom_hac_g303", 
        outdir=join(DATA_DIR, assembly_dir, "final_polished"),
    shell:
        """
        medaka_consensus -i {input.fq} -d {input.cn} -o {params.outdir} -t {threads} -m {params.model}
        mv {params.outdir}/consensus.fasta {output}
        """


rule pilon:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
        cn=join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{run}.medaka.fasta"),
        fw=ILLUMINA_FWD,
        rv=ILLUMINA_REV,
    output:
        join(DATA_DIR, assembly_dir, "final_polished/singlerun/{run}.polished.fasta"),
    threads: workflow.cores, 
    params:
        num_iter=2, 
        sample_name="{run}",
        outdir=join(DATA_DIR, assembly_dir, "final_polished/singlerun"),
    shell:
        """
        asm="{input.cn}"
        for i in {1..{params.num_iter}}
        do
            bwa index $asm
            bwa mem $asm {input.fw} {input.rv} | samtools view -bS - | samtools sort -@ {threads} - {params.outdir}/${i}.bam
            samtools index -@ {threads} {params.outdir}/${i}.bam
            pilon -Xmx$50g --threads {threads} --genome ${asm} --bam {params.outdir}/${i}.bam --output {params.outdir}/$i".pilon"
            asm={params.outdir}/$i".pilon.fasta"
        done
        mv {params.outdir}/"{params.num_iter}.pilon.fasta" {output}
        """


