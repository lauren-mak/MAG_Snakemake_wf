
from os.path import basename, exists
import pandas as pd


# From https://nf-co.re/mag/2.1.0/output 
rule porechop:
    input:
        join(FASTQ_DIR, "{sample}.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "porechop/{sample}.fastq.gz"),
    threads: workflow.cores,
    shell:
        """
        porechop -i {input} -t {threads} -o {output}
        """


def get_illumina_reads(sample, step):
    sample_file = "illumina.txt"
    df = pd.read_csv(sample_file, sep="\t")
    reads = df[df["sample"] == sample]["reads"].tolist()[0].split(",") # R1.fastq(.gz),R2.fastq(.gz)
    dict = {"fw": reads[0], "rv": reads[1]}
    print(sample)
    print(step)
    print(dict)
    return dict

# Parameters mostly from https://nf-co.re/mag/2.1.0/output 
# While nf-core MAG used Illumina reads to trim Nanopore reads, MUFFIN did not
rule filtlong:
    input:
        unpack(lambda wildcards: get_illumina_reads(wildcards.sample, "singlerun,filtlong")),
        fq=join(DATA_DIR, preprocessing_dir, "porechop/{sample}.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "filtlong/{sample}.fastq.gz"),
    params:
        short_qc="1000",
        long_len_wt="10",
    threads: workflow.cores,
    shell:
        """
        filtlong -1 {input.fw} -2 {input.rv} --trim  \
            --min_length {params.short_qc} \
            --mean_q_weight 10 \
            --split 100 \
            {input.fq} | gzip > {output}
        """


rule remove_host_genome:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "filtlong/{sample}.fastq.gz"),
        asm="/athena/masonlab/scratch/databases/metagenomics/HG37_Index/hg37dec_v0.1.fa",
    output:
        join(DATA_DIR, preprocessing_dir, "singlerun/{sample}.fastq.gz"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    params:
        sample_name="{sample}",
        outdir=join(DATA_DIR, preprocessing_dir, "singlerun/"),
    threads: workflow.cores,
    shell:
        """
        mkdir -p {params.outdir}
        prefix="{params.outdir}/{params.sample_name}.tmp"
        minimap2 -ax map-ont {input.asm} {input.fq} -o ${{prefix}}.sam
        samtools view -bS -f 4 ${{prefix}}.sam -o ${{prefix}}.bam
        samtools sort -@ {threads} -o ${{prefix}}.sort.bam ${{prefix}}.bam
        samtools fastq ${{prefix}}.sort.bam > ${{prefix}}.fastq
        gzip -c ${{prefix}}.fastq > {output}
        rm ${{prefix}}*
        """


# rule nanoplot:
#     input:
#         join(DATA_DIR, preprocessing_dir, "singlerun/{sample}.fastq.gz"),
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
        join(DATA_DIR, preprocessing_dir, "singlerun/{sample}.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "singlerun/{sample}_readcount.csv"),        
    params:
        sample_name="{sample}",
    shell:
        """
        num_reads=`cat {input} | awk 'NR%4==2{c++; l+=length($0)} END {print c","l}'`
        echo {sample_name},${num_reads} > {output}
        """
        # Adapted from https://wiki.itap.purdue.edu/display/CGSB/fastq%3A+count+reads+in+fastq+or+gzipped+fastq+file


rule count_reads_ss: # For the single-sample-only case
    input:
        join(DATA_DIR, binning_analyses, "singlerun/GTDB/gtdbtk.bac120.summary.tsv"), # Most effective way to ensure proper results-gathering
    output:
        join(DATA_DIR, preprocessing_dir, "readcounts.tsv"),        
    shell:
        """
        cat data/00_preprocessing/singlerun/*_readcount.csv > {output}
        """


def est_mem_flye(wildcards): # Unnecessary, for now
    info = open(join('data/00_preprocessing/singlerun', wildcards.sample + '_readcount.csv'),"r").readline().split(',')
    num_reads = int(info[2]) # sample_name,num_reads,num_bases
    base_mem = (int)((0.00000183542567428533 * num_reads - 8.01103718491264) * 1.1) # [-] Different?
    attempt_f = join(DATA_DIR, assembly_dir,'singlerun', wildcards.sample + '_attempt.txt') # Need this because wildcards.attempt is inaccessible?
    attempt_c = 0  
    if exists(attempt_f):
        attempt_c = (int)(open(attempt_f, 'r').readlines()[0].strip()) + 1
    print(attempt_c, file = open(attempt_f, 'w'))
    return max(max(base_mem, 0) + 100 * attempt_c, 50) # [-] Different?

rule metaflye:
    input:
        join(DATA_DIR, preprocessing_dir, "singlerun/{sample}.fastq.gz"),
        # join(DATA_DIR, preprocessing_dir, "singlerun/{sample}_readcount.txt"),        
    output:
        join(DATA_DIR, assembly_dir, "flye_unpolished/singlerun/{sample}/assembly.fasta"),
    threads: workflow.cores,
    params:
        outdir=join(DATA_DIR, assembly_dir, "flye_unpolished/singlerun/{sample}"),
    shell:
        """
        if [ -f "{params.outdir}/flye.log" && -f "{params.outdir}/00-assembly/draft_assembly.fasta" ]; then
            flye --resume -o {params.outdir} # Resume from the last previously completed step. Don't have to specify
        else
            flye --nano-raw {input} -o {params.outdir} -t {threads} --meta
        fi
        """


# Temporary BAM: {params.outdir}/{params.sample_name}_tmp.bam. Test to see if pipe works.
# rule polish_mapping:
#     input:
#         fq=join(DATA_DIR, preprocessing_dir, "singlerun/{sample}.fastq.gz"),
#         asm=join(DATA_DIR, assembly_dir, "flye_unpolished/singlerun/{sample}/assembly.fasta"),
#     output:
#         join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{sample}.bam"),
#     threads: workflow.cores,
#     params:
#         sample_name="{sample}",
#         outdir=join(DATA_DIR, assembly_dir, "prelim_polished/singlerun"),
#     shell:
#         """
#         minimap2 -ax map-ont {input.asm} {input.fq} | samtools view -bS - | samtools sort -@ {threads} -o {output} -
#         """

rule racon:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{sample}.fastq.gz"),
        asm=join(DATA_DIR, assembly_dir, "flye_unpolished/singlerun/{sample}/assembly.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{sample}/racon.fasta"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    params:
        num_iter=3, 
        outdir=join(DATA_DIR, assembly_dir, "prelim_polished/singlerun", "{sample}"),
    threads: workflow.cores,
    shell:
        """
        mkdir -p {params.outdir}
        asm="{input.asm}"
        for i in $( seq 1 "{params.num_iter}" )
        do
            prefix="{params.outdir}/${{i}}.tmp"
            minimap2 -ax map-ont $asm {input.fq} -o ${{prefix}}.sam
            samtools sort -@ {threads} -o ${{prefix}}.sort.sam ${{prefix}}.sam
            racon --t {threads} {input.fq} ${{prefix}}.sort.sam $asm > ${{prefix}}.racon.fasta
            asm="${{prefix}}.racon.fasta"
        done
        mv "{params.outdir}/{params.num_iter}.tmp.racon.fasta" {output}
        rm {params.outdir}/*tmp*
        """


rule medaka:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{sample}.fastq.gz"),
        cn=join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{sample}/racon.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{sample}/medaka.fasta"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    threads: workflow.cores,
    params:
        model="r941_prom_high_g4011", 
        outdir=join(DATA_DIR, assembly_dir, "prelim_polished/singlerun", "{sample}"),
    shell:
        """
        mkdir -p {params.outdir}
        medaka_consensus -i {input.fq} -d {input.cn} -o {params.outdir} -t {threads} -m {params.model}
        mv "{params.outdir}/consensus.fasta" {output}
        """


rule pilon:
    input:
        unpack(lambda wildcards: get_illumina_reads(wildcards.sample, "singlerun,pilon")),
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{sample}.fastq.gz"),
        cn=join(DATA_DIR, assembly_dir, "prelim_polished/singlerun/{sample}/medaka.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "final_polished/singlerun/{sample}.polished.fasta"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    threads: workflow.cores, 
    params:
        num_iter=2, 
        outdir=join(DATA_DIR, assembly_dir, "final_polished/singlerun", "{sample}"),
    shell:
        """
        mkdir -p {params.outdir}
        asm="{input.cn}"
        for i in $( seq 1 "{params.num_iter}" )
        do
            bwa index ${{asm}}
            prefix="{params.outdir}/${{i}}.tmp"
            bwa mem ${{asm}} {input.fw} {input.rv} -o ${{prefix}}.sam
            samtools view -bS -o ${{prefix}}.bam ${{prefix}}.sam 
            samtools sort -@ {threads} -o ${{prefix}}.sort.bam ${{prefix}}.bam
            samtools index -@ {threads} ${{prefix}}.sort.bam
            pilon -Xmx50g --threads {threads} --genome ${{asm}} --bam ${{prefix}}.sort.bam --output "${{prefix}}.pilon"
            asm="${{prefix}}.pilon.fasta"
        done
        mv "{params.outdir}/{params.num_iter}.tmp.pilon.fasta" {output}
        rm {params.outdir}/*tmp*
        """


