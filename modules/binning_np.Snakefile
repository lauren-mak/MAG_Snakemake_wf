
def get_sample_reads(wildcards):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    print(df, file=sys.stderr)
    datasets = df[df["coassembly"] == wildcards.sample]["datasets"][0].split(",") # S1,S2
    prefix = join(DATA_DIR, preprocessing_dir, "singlerun")
    return [join(prefix, i + ".fastq") for i in datasets]


def metawrap_cmmd(wildcards):
    return " ".join(get_sample_reads(wildcards))


checkpoint metawrap_binning:
    input:
        fastq=join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
        fasta=join(DATA_DIR, assembly_dir, "final_polished/singlerun/{run}.polished.fasta"),
    output:
        outfile=join(DATA_DIR, binning_dir, "singlerun/{run}/metawrap/metawrap_done.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/binning.yaml"
    params:
        fastq=join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq"),
        outdir=join(DATA_DIR, binning_dir, "singlerun/{run}/metawrap"),
        mincontiglength=1000,
    threads: workflow.cores
    resources:
        time=lambda wildcards, attempt: 20 * attempt,
        mem=lambda wildcards, attempt: 80 * attempt,
    shell:
        """
        rm -rf {params.outdir}
        gunzip {input.fastq} 
        metawrap binning -t {threads} -m {resources.mem} \
        -a {input.fasta} --maxbin2 --metabat2 \
        -l {params.mincontiglength} -o {params.outdir} --single-end {params.fastq}
        touch {output.outfile}
        """


checkpoint concoct_binning:
    input:
        join(DATA_DIR, binning_dir, "singlerun/{run}/metawrap/metawrap_done.txt"),
    output:
        outfile=join(DATA_DIR, binning_dir, "singlerun/{run}/metawrap/concoct_done.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/concoct.yaml" 
    params:
        outdir=join(DATA_DIR, binning_dir, "singlerun/{run}/metawrap"),
        mincontiglength=1000,
    threads: workflow.cores
    resources:
        time=lambda wildcards, attempt: 20 * attempt,
        mem=lambda wildcards, attempt: 80 * attempt,
    shell:
        """
        /home/lam4003/bin/MAG_Snakemake_wf/scripts/concoct.sh {params.outdir} {params.mincontiglength} 
        touch {output.outfile}
        """


checkpoint metawrap_binning_coas:
    input:
        unpack(get_sample_reads),
        join(DATA_DIR, assembly_dir, "final_polished/coassembly/{run}.polished.fasta"),
    output:
        outfile=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/metawrap_done.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/binning.yaml"
    params:
        assembly=join(DATA_DIR, assembly_dir, "final_polished/{run}.polished.fasta"),
        outdir=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/"),
        reads=lambda wildcards: metawrap_cmmd(wildcards),
        mincontiglength=1000,
    singularity:
        "shub://sskashaf/Containers:metawrap"
    threads: workflow.cores
    resources:
        time=lambda wildcards, attempt: 20 * attempt,
        mem=lambda wildcards, attempt: 80 * attempt,
    shell:
        """
        rm -rf {params.outdir}
        metawrap binning -t {threads} -m {resources.mem} \
        -a {params.assembly} --metabat2 --maxbin2  \
        -l {params.mincontiglength} -o {params.outdir} {params.reads}
        touch {output}
        """


checkpoint concoct_binning_coas:
    input:
        join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/metawrap_done.txt"),
    output:
        outfile=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/concoct_done.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/concoct.yaml"
    params:
        outdir=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/"),
        mincontiglength=1000,
    threads: workflow.cores
    resources:
        time=lambda wildcards, attempt: 20 * attempt,
        mem=lambda wildcards, attempt: 80 * attempt,
    shell:
        """
        /home/lam4003/bin/MAG_Snakemake_wf/scripts/concoct.sh {params.outdir} {params.mincontiglength}
        touch {output.outfile}
        """


