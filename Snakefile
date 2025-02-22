# -*- coding: utf-8

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

#'''
#   This is a basic framework for recovery and basic quality control of MAGs
#   To visualize the pipeline: snakemake --dag | dot -Tpng > dag.png
#'''

__maintainer__ = "Sara Kashaf"
__email__ = "sskashaf@ebi.ac.uk"


import os
from os.path import join
import sys
import glob
import pandas as pd
import csv

configfile: "config.yaml"

# Directory structure
DATA_DIR = "data" #config['data']
preprocessing_dir = "00_preprocessing"
assembly_dir = "01_assembly"
binning_dir = "02_binning"
binning_analyses = "03_binning_analyses"
if not os.path.exists("logs"):
    os.makedirs("logs")
os.system("chmod +x scripts/*")
os.system("chmod +x scripts/plotting/*")

# LOAD METADATA
df_run = pd.read_csv("runs.txt")
RUN = df_run["Run"]

df_coas = pd.read_csv("coassembly_runs.txt", sep="\t")
COAS = df_coas["coassembly"]


all_outfiles = [
     # Figure 2
     join(DATA_DIR, "final_reports", "raw_multiqc_report.html"),
     join(DATA_DIR, "final_reports", "post_multiqc_report.html"),
     # Figure 3a
     join(DATA_DIR, "final_reports", "sr_summ_cmseq_all.csv"),
     join(DATA_DIR, "final_reports", "coas_summ_cmseq_all.csv"),
     join(DATA_DIR, "final_reports", "sr_checkm_metrics.csv"),
     join(DATA_DIR, "final_reports", "coas_checkm_metrics.csv"),
     # Figure 3b
     join(DATA_DIR, "final_reports", "dnadiff_summary.tsv"),
     # Figure 3c
     join(DATA_DIR, "final_reports", "gtdbtk.bac120.summary.tsv"),
     # Figure 4
     join(DATA_DIR, "final_reports", "readcounts.tsv"),
     join(DATA_DIR, "final_reports", "flagstat_sum.txt"),
     join(DATA_DIR, "final_reports", "sr_catalogue_mapreads.tab"),
     join(DATA_DIR, "final_reports", "coas_catalogue_mapreads.tab"),
]

rule copy_output:
    input:
        rawd_qc_in=join(DATA_DIR,preprocessing_dir, "raw_qc/multiqc/raw_multiqc_report.html"),
        post_qc_in=join(DATA_DIR,preprocessing_dir, "postprocessing_qc/multiqc/post_multiqc_report.html"),
        cmseq_1_in=join(DATA_DIR, binning_analyses, "singlerun/cmseq/summ_cmseq_all.csv"),
        cmseq_2_in=join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/summ_cmseq_all.csv"),
        gtdbt_1_in=join(DATA_DIR, binning_analyses, "singlerun_coassembly/GTDB/gtdbtk.bac120.summary.tsv"),
        fwork_1_in=join(DATA_DIR, preprocessing_dir, "readcounts.tsv"),
        fwork_2_in=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/flagstat_sum.txt"),
        fwork_3_in=join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/sr_catalogue_mapreads.tab"),
        fwork_4_in=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/coas_catalogue_mapreads.tab"),
        chckm_1_in=join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics.csv"),
        chckm_2_in=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.csv"),
        ddiff_3_in=join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff/dnadiff_summary.tsv"),
    output:
        rawd_qc_out=join(DATA_DIR, "final_reports", "raw_multiqc_report.html"),
        post_qc_out=join(DATA_DIR, "final_reports", "post_multiqc_report.html"),
        cmseq_1_out=join(DATA_DIR, "final_reports", "sr_summ_cmseq_all.csv"),
        cmseq_2_out=join(DATA_DIR, "final_reports", "coas_summ_cmseq_all.csv"),
        gtdbt_1_out=join(DATA_DIR, "final_reports", "gtdbtk.bac120.summary.tsv"),
        fwork_1_out=join(DATA_DIR, "final_reports", "readcounts.tsv"),
        fwork_2_out=join(DATA_DIR, "final_reports", "flagstat_sum.txt"),
        fwork_3_out=join(DATA_DIR, "final_reports", "sr_catalogue_mapreads.tab"),
        fwork_4_out=join(DATA_DIR, "final_reports", "coas_catalogue_mapreads.tab"),
        chckm_1_out=join(DATA_DIR, "final_reports", "sr_checkm_metrics.csv"),
        chckm_2_out=join(DATA_DIR, "final_reports", "coas_checkm_metrics.csv"),
        ddiff_3_out=join(DATA_DIR, "final_reports", "dnadiff_summary.tsv"),
    params:
        outdir=join(DATA_DIR, "final_reports"),
    shell:
        """
        mkdir -p {params.outdir}
        cp {input.rawd_qc_in} {output.rawd_qc_out} 
        cp {input.post_qc_in} {output.post_qc_out} 
        cp {input.cmseq_1_in} {output.cmseq_1_out} 
        cp {input.cmseq_2_in} {output.cmseq_2_out} 
        cp {input.gtdbt_1_in} {output.gtdbt_1_out} 
        cp {input.fwork_1_in} {output.fwork_1_out} 
        cp {input.fwork_2_in} {output.fwork_2_out} 
        cp {input.fwork_3_in} {output.fwork_3_out} 
        cp {input.fwork_4_in} {output.fwork_4_out} 
        cp {input.chckm_1_in} {output.chckm_1_out} 
        cp {input.chckm_2_in} {output.chckm_2_out} 
        cp {input.ddiff_3_in} {output.ddiff_3_out} 
        """

rule all:
   input: all_outfiles


include: "modules/sra_download.Snakefile"
include: "modules/preprocessing.Snakefile"
include: "modules/coas.Snakefile"
include: "modules/assembly.Snakefile"
include: "modules/binning.Snakefile"
include: "modules/refine.Snakefile"
include: "modules/dRep_GTDB.Snakefile"
include: "modules/framework.Snakefile"
include: "modules/refine_coas.Snakefile"
include: "modules/cmseq.Snakefile"
include: "modules/dnadiff.Snakefile"
