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

# vim: set ft=python

#'''
#   Mod purpose: Ingests local FastQ files instead of downloads from SRA prefetch. Applies BayesHammer before SPAdes step.
#   Author: Lauren Mak
#   Last updated: 10/6/2021
#'''

from os.path import basename
import pandas as pd

def est_mem_spades(wildcards):
    read_counts = pd.read_csv('data/00_preprocessing/readcounts.tsv', sep = '\t', header = 0)
    read_counts['Run'] = ['_'.join(basename(i).split('_')[0:3]) for i in read_counts['Run']]

    singlerun_reads = read_counts[read_counts['Run'] == wildcards.run]
    num_reads = singlerun_reads['Readcount'].sum()
    base_mem = (int)((0.00000183542567428533 * num_reads - 8.01103718491264) * 1.1)
    attempt_f = join(DATA_DIR, assembly_dir,'singlerun', wildcards.run + '_attempt.txt') # Need this because wildcards.attempt is inaccessible?
    attempt_c = 0  
    if exists(attempt_f):
        attempt_c = (int)(open(attempt_f, 'r').readlines()[0].strip()) + 1
    print(attempt_c, file = open(attempt_f, 'w'))
    return max(max(base_mem, 0) + 100 * attempt_c, 50) * 1000

checkpoint spades:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_2.fastq.gz"),
        upd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_upd.fastq.gz"),
    output:
        join(DATA_DIR, assembly_dir, "singlerun/{run}/scaffolds.fasta"),
        join(DATA_DIR, assembly_dir, "singlerun/{run}/contigs.fasta"),
    threads: workflow.cores,
    params:
        outdir=join(DATA_DIR, assembly_dir, "singlerun/{run}/"),
    # resources:
    #     time=lambda wildcards, attempt: 10 * attempt,
    #     mem=est_mem_spades, # lambda wildcards, attempt: 200 + 100 * attempt,
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:assembly"
    shell:
        """
        # rm -rf {params.outdir}
        # dmesg -T
        if [ -f "{params.outdir}/K33/assembly_graph.fastg" ]; then
            spades.py --restart-from k55 -o {params.outdir} -t {threads} -m {resources.mem} # Restart from k = 55
        elif [ -f "{params.outdir}/K21/assembly_graph.fastg" ]; then
            spades.py --restart-from k33 -o {params.outdir} -t {threads} -m {resources.mem} # Restart from k = 33
        elif [ -f "{params.outdir}/corrected.yaml" ]; then
            spades.py --restart-from k21 -o {params.outdir} -t {threads} -m {resources.mem} # Restart from k = 21
        else
            spades.py --only-assembler --meta -1 {input.fwd} -2 {input.rev} -s {input.upd} -o {params.outdir} -t {threads} -m {resources.mem}
        fi
        """

def est_mem_coas(wildcards):
    read_counts = pd.read_csv('data/00_preprocessing/readcounts.tsv', sep = '\t', header = 0)

    num_reads = read_counts['Readcount'].sum()
    base_mem = (int)((0.00000183542567428533 * num_reads - 8.01103718491264) * 1.1)
    attempt_f = join(DATA_DIR, assembly_dir,'coassembly', wildcards.run + '_attempt.txt') # Need this because wildcards.attempt is inaccessible?
    attempt_c = 0  
    if exists(attempt_f):
        attempt_c = (int)(open(attempt_f, 'r').readlines()[0].strip()) + 1
    print(attempt_c, file = open(attempt_f, 'w'))
    return max(max(base_mem, 0) + 50 * attempt_c, 50) * 1000

checkpoint spades_coas:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_2.fastq.gz"),
        upd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_upd.fastq.gz"),
    output:
        join(DATA_DIR, assembly_dir, "coassembly/{run}/scaffolds.fasta"),
        join(DATA_DIR, assembly_dir, "coassembly/{run}/contigs.fasta"),
    # resources:
    #     time=lambda wildcards, attempt: 10 * attempt,
    #     mem=est_mem_coas, # lambda wildcards, attempt: 200 + 150 * attempt,
    threads: workflow.cores,
    params:
        outdir=join(DATA_DIR, assembly_dir, "coassembly/{run}/"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:assembly"
    shell:
        """
        # rm -rf {params.outdir}
        # dmesg -T
        if [ -f "{params.outdir}/K33/assembly_graph.fastg" ]; then
            spades.py --restart-from k55 -o {params.outdir} -t {threads} -m {resources.mem} # Restart from k = 55
        elif [ -f "{params.outdir}/K21/assembly_graph.fastg" ]; then
            spades.py --restart-from k33 -o {params.outdir} -t {threads} -m {resources.mem} # Restart from k = 33
        elif [ -f "{params.outdir}/corrected.yaml" ]; then
            spades.py --restart-from k21 -o {params.outdir} -t {threads} -m {resources.mem} # Restart from k = 21
        else
            spades.py --only-assembler --meta -1 {input.fwd} -2 {input.rev} -s {input.upd} -o {params.outdir} -t {threads} -m {resources.mem}
        fi
        """