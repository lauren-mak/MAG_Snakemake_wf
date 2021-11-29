#!/bin/bash

# $1 = sample name

# Set FastQ variables
fw_fastq="data/00_preprocessing/kneaddata_bowtie/singlerun/${1}_1.fastq"
rv_fastq="data/00_preprocessing/kneaddata_bowtie/singlerun/${1}_2.fastq"
interleaved_fastq="data/00_preprocessing/processed/${1}/interleaved.fastq"
processed_fastq="data/00_preprocessing/processed/${1}/processed.fastq"
corrected_fastq="data/00_preprocessing/processed/${1}/corrected.fastq"

# Interleave FastQs. Adapted from https://gist.github.com/nathanhaigh/4544979
# Must preprocess anyways because there are still Ns in the reads and index can't handle that
# Must be done before loading GLIBC! If not, Java will hang
mkdir -p "data/00_preprocessing/processed/${sample}"
paste ${fw_fastq} ${rv_fastq} | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}' > ${interleaved_fastq} &

export LD_LIBRARY_PATH=/home/lam4003/bin/anaconda3/envs/snakemake/glibc-2.14/lib

# Eliminate Ns and interleave FastQs 
sga preprocess --permute-ambiguous --no-primer-check --pe-mode=2 --pe-orphans="data/00_preprocessing/processed/singlerun/${sample}_upd.fastq" --out=${processed_fastq} ${interleaved_fastq}

# Index temporary interleaved FastQ 
sga index -a ropebwt -t 30 --no-reverse ${processed_fastq}

# Correct interleaved FastQ 
sga correct --learn -t 30 --out=${corrected_fastq} ${processed_fastq}

# Deinterleave FastQ. Adapted from https://gist.github.com/nathanhaigh/3521724
cat ${corrected_fastq} | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | \
    gzip --best > "data/00_preprocessing/processed/singlerun/${sample}_1.fastq.gz") | cut -f 5-8 | tr "\t" "\n" | \
    gzip --best > "data/00_preprocessing/processed/singlerun/${sample}_2.fastq.gz"

# Remove intermediate files 
rm ${interleaved_fastq} ${processed_fastq} ${corrected_fastq}
# Make an unpaired FastQ to match assembly input
touch "data/00_preprocessing/processed/singlerun/${sample}_upd.fastq"
gzip "data/00_preprocessing/processed/singlerun/${sample}_upd.fastq"
