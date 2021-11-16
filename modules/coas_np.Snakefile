
def get_sample_reads_coas(sample):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    print(df, file=sys.stderr)
    datasets = df[df["coassembly"] == sample]["datasets"][0].split(",") # S1,S2
    prefix = join(DATA_DIR, preprocessing_dir, "singlerun")
    return [join(prefix, i + ".fastq.gz") for i in datasets]


rule merge_coas:
    input:
        unpack(lambda wildcards: get_sample_reads_coas(wildcards.sample)),
    output:
        join(DATA_DIR, preprocessing_dir, "coassembly/{sample}.fastq.gz"),
    shell:
        """
        cat {input} > {output}
        """


rule count_reads_coas:
    input:
        join(DATA_DIR, preprocessing_dir, "coassembly/{sample}.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "readcounts.tsv"),        
    shell:
        """
        cat data/00_preprocessing/singlerun/*_readcount.csv > {output}
        """


rule metaflye_coas:
    input:
        join(DATA_DIR, preprocessing_dir, "coassembly/{sample}.fastq.gz"),
    output:
        join(DATA_DIR, assembly_dir, "flye_unpolished/coassembly/{sample}/assembly.fasta"),
    threads: workflow.cores,
    params:
        outdir=join(DATA_DIR, assembly_dir, "flye_unpolished/coassembly/{sample}"),
    shell:
        """
        if [ -f "{params.outdir}/flye.log" ]; then
            flye --resume -o {params.outdir} # Resume from the last previously completed step. Don't have to specify
        else
            flye --nano-raw {input} -o {params.outdir} -t {threads} --meta
        fi
        """

# Temporary BAM: {params.outdir}/{params.sample_name}_tmp.bam. Test to see if pipe works.
rule polish_mapping_coas:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "coassembly/{sample}.fastq.gz"),
        asm=join(DATA_DIR, assembly_dir, "flye_unpolished/coassembly/{sample}/assembly.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "prelim_polished/coassembly/{sample}.bam"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    threads: workflow.cores,
    params:
        sample_name="{sample}",
        outdir=join(DATA_DIR, assembly_dir, "prelim_polished/coassembly"),
    shell:
        """
        minimap2 -ax map-ont {input.asm} {input.fq} | samtools view -bS - | samtools sort -@ {threads} -o {output} -
        """

rule racon_coas:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "coassembly/{sample}.fastq.gz"),
        mmp=join(DATA_DIR, assembly_dir, "prelim_polished/coassembly/{sample}.bam"),
        asm=join(DATA_DIR, assembly_dir, "flye_unpolished/coassembly/{sample}/assembly.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "prelim_polished/coassembly/{sample}.racon.fasta"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    threads: workflow.cores,
    params:
        num_iter=3, 
        outdir=join(DATA_DIR, assembly_dir, "prelim_polished/coassembly"),
    threads: workflow.cores,
    shell:
        """
        asm="{input.asm}"
        for i in {1..{params.num_iter}}
        do
            minimap2 -ax map-ont $asm {input.fq} | samtools view -bS - | samtools sort -@ {threads} -o {params.outdir}/${i}.bam -
            samtools index -@ {threads} {params.outdir}/${i}.bam
            racon --t {threads} {input.fq} {params.outdir}/${i}.bam $asm > {params.outdir}/${i}.racon.fasta
            asm={params.outdir}/${i}.racon.fasta
        done
        mv {params.outdir}/"{params.num_iter}.racon.fasta" {output}
        """


rule medaka_coas:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "coassembly/{sample}.fastq.gz"),
        cn=join(DATA_DIR, assembly_dir, "prelim_polished/coassembly/{sample}.racon.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "prelim_polished/coassembly/{sample}.medaka.fasta"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    threads: workflow.cores,
    params:
        model="r941_prom_hac_g303", 
        outdir=join(DATA_DIR, assembly_dir, "final_polished"),
    shell:
        """
        medaka_consensus -i {input.fq} -d {input.cn} -o {params.outdir} -t {threads} -m {params.model}
        mv {params.outdir}/consensus.fasta {output}
        """


def get_illumina_reads_coas(sample):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    print(df, file=sys.stderr)
    datasets = df[df["coassembly"] == sample]["datasets"][0].split(",") # S1,S2
    fw = []
    rv = []
    for i in dataset:
        r = get_illumina_reads(i)
        fw.append(r["fw"]) # [x_1.fastq, x_2.fastq, ...]
        rv.append(r["rv"])
    dict = {"fw": " ".join(fw), "rv": " ".join(rv)} # {fw: "x_1.fastq x_2.fastq"}
    print(dict)
    return dict


rule pilon_coas:
    input:
        unpack(lambda wildcards: get_illumina_reads(wildcards.sample)),
        fq=join(DATA_DIR, preprocessing_dir, "coassembly/{sample}.fastq.gz"),
        cn=join(DATA_DIR, assembly_dir, "prelim_polished/coassembly/{sample}.medaka.fasta"),
    output:
        join(DATA_DIR, assembly_dir, "final_polished/coassembly/{sample}.polished.fasta"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    threads: workflow.cores, 
    params:
        num_iter=2, 
        sample_name="{sample}",
        outdir=join(DATA_DIR, assembly_dir, "final_polished/coassembly"),
    shell:
        """
        asm="{input.cn}"
        mkdir -p {params.outdir}
        cat {input.fw} > {params.outdir}/tmp_1.fastq
        cat {input.rv} > {params.outdir}/tmp_2.fastq
        for i in {1..{params.num_iter}}
        do
            bwa index $asm
            bwa mem $asm {params.outdir}/tmp_1.fastq {params.outdir}/tmp_2.fastq | samtools view -bS - | samtools sort -@ {threads} - {params.outdir}/${i}.bam
            samtools index -@ {threads} {params.outdir}/${i}.bam
            pilon -Xmx$50g --threads {threads} --genome ${asm} --bam {params.outdir}/${i}.bam --output {params.outdir}/$i".pilon"
            asm={params.outdir}/$i".pilon.fasta"
        done
        mv {params.outdir}/"{params.num_iter}.pilon.fasta" {output}
        rm {params.outdir}/tmp*
        """


