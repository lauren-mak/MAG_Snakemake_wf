# This file is part of MAG Snakemake workflow.
#
# MAG Snakemake workflow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MAG Snakemake workflow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with MAG Snakemake workflow.  If not, see <https://www.gnu.org/licenses/>.

# vim::w set ft=python:

from os.path import basename, exists
import pandas as pd

rule raw_fastqc_fwd:
    input:
        fwd=join(FASTQ_DIR, "{run}_1.fastq.gz"), # New
        rev=join(FASTQ_DIR, "{run}_2.fastq.gz"), # New
    output:
        join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/{run}_1_fastqc.html"),
    params:
        outdir=directory(join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/")),
    threads: workflow.cores
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    shell:
        """
        fastqc {input.fwd} --outdir {params.outdir}
        """


rule raw_fastqc_rev:
    input:
        fwd=join(FASTQ_DIR, "{run}_1.fastq.gz"), # New
        rev=join(FASTQ_DIR, "{run}_2.fastq.gz"), # New
    output:
        join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/{run}_2_fastqc.html"),
    params:
        outdir=directory(join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/")),
    threads: workflow.cores
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    shell:
        """
        fastqc {input.rev} --outdir {params.outdir}
        """



rule raw_multiqc:
    input:
        expand(join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/{run}_{read}_fastqc.html"), run=RUN, read=["1", "2"]),
    output:
        join(DATA_DIR, preprocessing_dir, "raw_qc/multiqc/raw_multiqc_report.html"),
    params:
        indir=join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/"),
        outdir=join(DATA_DIR, preprocessing_dir, "raw_qc/multiqc/"),
        outfile=join(DATA_DIR, preprocessing_dir, "raw_qc/multiqc/multiqc_report.html"),
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.3--py35_2"
    shell:
        """
        rm -f {params.outfile}
        multiqc {params.indir} -o {params.outdir}
        mv {params.outfile} {output}
        """


rule kneaddata_download_database:
    output:
        join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.1.bt2"),
        join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.2.bt2"),
        join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.3.bt2"),
        join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.4.bt2"),
        join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.rev.1.bt2"),
        join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.rev.2.bt2"),
    params:
        outdir="/athena/masonlab/scratch/databases/metagenomics/HG37_Index/",
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:metagenome_preprocessing"
    shell:
        """
        kneaddata_database --download human_genome bowtie2 {params.outdir}
        """

rule kneaddata_bowtie:
    input:
        fwd=join(FASTQ_DIR, "{run}_1.fastq.gz"), # New
        rev=join(FASTQ_DIR, "{run}_2.fastq.gz"), # New
        indx1=join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.1.bt2"),
        indx2=join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.2.bt2"),
        indx3=join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.3.bt2"),
        indx4=join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.4.bt2"),
        indx5=join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.rev.1.bt2"),
        indx6=join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/hg37dec_v0.1.rev.2.bt2"),
    output:
        fwd=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_2.fastq"),
    params:
        tmp_fwd=join(DATA_DIR, "{run}_tmp_1.fastq"),
        tmp_rev=join(DATA_DIR, "{run}_tmp_2.fastq"),
        tmp_fwd2=join(DATA_DIR, "{run}_tmp2_1.fastq.gz"),
        tmp_rev2=join(DATA_DIR, "{run}_tmp2_2.fastq.gz"),
        fwd=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/{run}_tmp2_1_kneaddata_paired_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/{run}_tmp2_1_kneaddata_paired_2.fastq"),
        outdir=directory(join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/")),
        indx=join("/athena/masonlab/scratch/databases/metagenomics", "HG37_Index/"),
        prddir=directory(join(DATA_DIR, preprocessing_dir, "processed/")),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:metagenome_preprocessing"
    threads: workflow.cores,
    shell:
        """
        #this version of kneaddata requires read identifiers (/1, /2). Given we did not download the file from the sra using the option --readids, we need to remove spaces in the headers
        # and then add /1 and /2 using reformat.sh. Kneaddata by default does not sort the reads so Next we sort the reads. 
        # Alternatively, one could sort using kneaddata with the option --reorder
        rm -f {params.tmp_fwd} {params.tmp_rev} {params.tmp_fwd2} {params.tmp_rev2} {params.fwd} {params.rev} {output.fwd} {output.rev}
        mkdir -p {params.outdir}
        mkdir -p {params.prddir}
        seqtk seq -C {input.fwd} > {params.tmp_fwd}
        seqtk seq -C {input.rev} > {params.tmp_rev}
        /home/lam4003/bin/MAG_Snakemake_wf/scripts/reformat.sh in={params.tmp_fwd} in2={params.tmp_rev} out1={params.tmp_fwd2} out2={params.tmp_rev2} addslash spaceslash=f
        kneaddata --remove-intermediate-output --threads {threads} \
        --input {params.tmp_fwd2} --input {params.tmp_rev2}\
        --output {params.outdir} \
        --reference-db {params.indx} \
        --trimmomatic-options "ILLUMINACLIP:/data/adapters/TruSeq3-PE.fa:2:30:10: SLIDINGWINDOW:4:20 MINLEN:50" --trimmomatic /home/lam4003/bin/anaconda3/share/trimmomatic-0.39-2\
        --bowtie2-options "--very-sensitive --dovetail" 
        /home/lam4003/bin/bbmap/repair.sh in={params.fwd} in2={params.rev} out={output.fwd} out2={output.rev} repair
        rm {params.tmp_fwd} {params.tmp_rev} {params.tmp_fwd2} {params.tmp_rev2}
        """


def est_mem_spades(wildcards):
    read_counts = pd.read_csv('data/00_preprocessing/readcounts.tsv', sep = '\t', header = 0)
    read_counts['Run'] = ['_'.join(basename(i).split('_')[:-1]) for i in read_counts['Run']]

    singlerun_reads = read_counts[read_counts['Run'] == wildcards.run]
    num_reads = singlerun_reads['Readcount'].sum()
    base_mem = (int)((0.00000183542567428533 * num_reads - 8.01103718491264) * 1.1)
    attempt_f = join(DATA_DIR, preprocessing_dir,'processed', wildcards.run + '_attempt.txt') # Need this because wildcards.attempt is inaccessible?
    attempt_c = 0  
    if exists(attempt_f):
        attempt_c = (int)(open(attempt_f, 'r').readlines()[0].strip()) + 1
    print(attempt_c, file = open(attempt_f, 'w'))
    return max(max(base_mem, 0) + 100 * attempt_c, 50) * 1000


rule bayeshammer:
    input:
        rct=join(DATA_DIR, preprocessing_dir, "readcounts.tsv"),
        fwd=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_2.fastq"),
    output:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_2.fastq.gz"),
        upd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_upd.fastq.gz"),
    params:
        prefix=join(DATA_DIR, preprocessing_dir, "processed", "{run}", "corrected", "{run}"),
        outdir=directory(join(DATA_DIR, preprocessing_dir, "processed", "{run}")),
    threads: workflow.cores,
    # resources:
    #     mem=est_mem_spades,
    shell:
        """
        mkdir -p data/01_assembly/singlerun
        mkdir -p data/01_assembly/coassembly
        spades.py --only-error-correction --meta -1 {input.fwd} -2 {input.rev} -o {params.outdir} -t {threads} -m 1000
        mv {params.prefix}_1.00.0_0.cor.fastq.gz {output.fwd}
        mv {params.prefix}_2.00.0_0.cor.fastq.gz {output.rev}
        mv {params.prefix}__unpaired.00.0_0.cor.fastq.gz {output.upd}
        """


rule postpreprocessing_fastqc_fwd:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_2.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/{run}_1_fastqc.html"),
    params:
        outdir=directory(join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/")),
    threads: workflow.cores,
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    shell:
        """
        fastqc {input.fwd} --outdir {params.outdir}
        """

rule postpreprocessing_fastqc_rev:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_2.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/{run}_2_fastqc.html"),
    params:
        outdir=directory(join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/")),
    threads: workflow.cores
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    shell:
        """
        fastqc {input.rev} --outdir {params.outdir}
        """



rule postpreprocessing_multiqc:
    input:
        expand(join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/{run}_{read}_fastqc.html"), run=RUN, read=["1", "2"]),
    output:
        join(DATA_DIR, preprocessing_dir, "postprocessing_qc/multiqc/post_multiqc_report.html"),
    params:
        indir=join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/"),
        outdir=join(DATA_DIR, preprocessing_dir, "postprocessing_qc/multiqc/"),
        outfile=join(DATA_DIR, preprocessing_dir, "postprocessing_qc/multiqc/multiqc_report.html"),
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.3--py35_2"
    shell:
        """
        multiqc --force {params.indir} -o {params.outdir}
        mv {params.outfile} {output}
        """


def linecount(fastq):
    out = subprocess.Popen("zcat -f "+fastq+" | wc -l", stdout=subprocess.PIPE, shell=True).communicate()[0]
    count=int(out.strip())/4
    return count

rule readcount_fq:
    input:
        raw=expand(join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_{read}.fastq"), run=RUN, read=["1", "2"]),
    output:
        join(DATA_DIR, preprocessing_dir, "readcounts.tsv"),
    run:
        outfile = str(output)
        if os.path.exists(outfile):
            os.remove(outfile)
        with open(outfile, "w") as outf:
            outf.writelines("Run\tReadcount\n")
            for run in input:
                readcount = int(linecount(run))
                line = run + "\t" + str(readcount)
                outf.writelines(line + "\n")


