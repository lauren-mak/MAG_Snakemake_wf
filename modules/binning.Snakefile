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

# vim: set ft=python:


def get_sample_reads_1(wildcards):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    print(df, file=sys.stderr)
    datasets = df[df["coassembly"] == wildcards.sample]["datasets"][0].split(",") # S1,S2
    prefix = join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie", "singlerun")
    r1_lst = [join(prefix, i + "_1.fastq") for i in datasets]
    return r1_lst


def get_sample_reads_2(wildcards):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    print(df, file=sys.stderr)
    datasets = df[df["coassembly"] == wildcards.sample]["datasets"][0].split(",") # S1,S2
    prefix = join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie", "singlerun")
    r2_lst = [join(prefix, i + "_2.fastq") for i in datasets]
    return r2_lst


def metawrap_cmmd(wildcards):
    r1_lst = get_sample_reads_1(wildcards)
    r2_lst = get_sample_reads_2(wildcards)
    cmmd = " ".join(r1_lst) + " " + " ".join(r2_lst)
    return cmmd


checkpoint metawrap_binning:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{sample}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{sample}_2.fastq"),
        fasta=join(DATA_DIR, assembly_dir, "singlerun/{sample}/scaffolds.fasta"),
    output:
        outfile=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap/metawrap_done.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/binning.yaml"
    singularity:
        "shub://sskashaf/Containers:metawrap"
    params:
        outdir=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap"),
        mincontiglength=1000,
    threads: workflow.cores
    resources:
        time=lambda wildcards, attempt: 20 * attempt,
        mem=lambda wildcards, attempt: 80 * attempt,
    shell:
        """
        rm -rf {params.outdir}
        metawrap binning -t {threads} -m {resources.mem}\
        -a {input.fasta} --maxbin2 --metabat2  \
        -l {params.mincontiglength} -o {params.outdir} {input.fwd} {input.rev}
        touch {output.outfile}
        """


checkpoint concoct_binning:
    input:
        join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap/metawrap_done.txt"),
    output:
        outfile=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap/concoct_done.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/concoct.yaml" 
    params:
        outdir=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap"),
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
        unpack(get_sample_reads_1),
        unpack(get_sample_reads_2),
        join(DATA_DIR, assembly_dir, "coassembly/{sample}/scaffolds.fasta"),
    output:
        outfile=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/metawrap_done.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/binning.yaml"
    params:
        assembly=join(DATA_DIR, assembly_dir, "coassembly/{sample}/scaffolds.fasta"),
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
        mincontiglength=2500,
    threads: workflow.cores
    resources:
        time=lambda wildcards, attempt: 20 * attempt,
        mem=lambda wildcards, attempt: 120 * attempt,
    shell:
        """
        /home/lam4003/bin/MAG_Snakemake_wf/scripts/concoct.sh {params.outdir} {params.mincontiglength}
        touch {output.outfile}
        """


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

