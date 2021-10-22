# MAG Snakemake Workflow

## Contents

- **Pipeline Setup**: Installing the pipeline and the required Snakemake conda environment.
- **Running the Pipeline**: Requirements for input files and description of output files. Sample submission commands and scripts for both local and scheduled jobs.
- **Modifications**: Summary of modifications that have been made to the original pipeline.
- **More Information**

## Pipeline Setup

### MAG Snakemake Workflow

```Bash
git clone https://github.com/lauren-mak/MAG_Snakemake_wf.git
```
- An Athena-optimized version of Saheb Kashaf et al. 2021's MAG Snakemake Workflow

### Conda Environments

- I found that four were required to separate packages with conflicting dependencies. 
    - Snakemake: The main environment, houses most of the programs.
    - Binning: The only unbroked conda version of MetaWRAP relies on Python 2, which is incompatible with Snakemake.
    - CONCOCT: GLIBC_2.14 issues with the most recent version, and environment incompatibilities with older versions. Easiest way around it was just to silo it in its own environment.
    - Prokka: Only for running Prokka. Older versions depended on a broken libcurl library from conda, and removing Prokka from the `snakemake` and `binning` environments was too difficult. 
- These are already set up at `/home/lam4003/bin/anaconda3/envs` and the Snakefile rules directly refer to these, so you don't have to set them up yourself if you're also working in Athena, **with the exception of the main Snakemake environment**. The installation instructions for that are below. 

#### Snakemake

```Bash
conda create -c conda-forge bioconda -n snakemake snakemake=5.18
conda activate snakemake
conda install -c bioconda parallel-fastq-dump=0.6.7 fastqc multiqc kneaddata seqtk spades crimson checkm-genome cmseq mummer drep=3.2.2 gtdbtk bwa samtools=1.12 mash
```
### Databases

- Location: `/athena/masonlab/scratch/databases/metagenomics`

Database | Recommended Version | Downloaded Version | Notes
-------- | ------------------- | ------------------ | -----
RefSeq complete bacterial genomes | May 2020         | January 8, 2020 (Release 98) | NA
GTDB database                     | Release 95       | Release 202      | NA
CheckM                            | January 16, 2015 | January 16, 2015 | This is the only version available

- **These are already set up if you're also working in Athena**. 
- For more information on how each of these databases were downloaded and the RefSeq Mash sketch made, DM me directly.

## Running the Pipeline

### Input

- The following files must be copied from `~/bin/MAG_Snakemake_wf/` to the working directory (ex. `working_dir/`):
- `runs.txt`: Each row contains a sample prefix. Corresponds to `X` as described in 'Input FastQs'.
    - If direct-downloading from the SRA, the sample prefix is the SRA accession ID of each sample.
- `coassembly_runs.txt`: Must be a **TSV** file. First column is the coassembled set name. Second column is a comma-separated list of sample prefixes.
- `config.yaml`: Paths to the GTDB and RefSeq databases.
- `clusterconfig.yaml`: Slurm submission settings. 
- `Snakefile*`: Set the variable `FASTQ_DIR` to the base location of the input FastQs.
    - `Snakefile`: For FastQs directly downloaded from the SRA.
    - `Snakefile_LD_SC`: For i) FastQs available locally to be ii) single-sample and co-assembled. Must be renamed to 'Snakefile'.
    - `Snakefile_LD_SS`: For i) FastQs available locally to be ii) single-sample-assembled only. Must be renamed to 'Snakefile'. For coassembly-only, `cat` all single-sample FastQs together before using the workflow before treating as a single-sample-only run.

#### Input FastQs 

- Gzipped and located in the same directory (does not need to be the working directory)
- Format: `X_N.fastq.gz` , where `X` is the sample prefix and `N` is the paired-end direction (1 or 2). 
    - Underscores (`_`) may be used in the sample prefix `X`

### Output

### Final Reports

- All output and intermediate files will be found in `working_dir/data`. If the job was Slurm-scheduled, the rule-speciic logs can be found in `working_dir/logs`. If the job was run locally, no rule-specific logs, unfortunately.
- The following reports will be found in `working_dir/data/final_reports`:
    - **raw_multiqc_report.html**: MultiQC report of all input FastQs before Kneaddata cleaning.
    - **post_multiqc_report.html**: MultiQC report of all input FastQs after Kneaddata cleaning.
    - **sr_summ_cmseq_all.csv**: CMSeq (strain heterogeneity) summary statistics of single-sample MAGs. 
    - **coas_summ_cmseq_all.csv**: CMSeq (strain heterogeneity) summary statistics of coassembled MAGs. 
    - **gtdbtk.bac120.summary.tsv**: GTDB classifications of each classifiable refined MAG bin, with the ANI to the most similar RefSeq reference genome. 
    - **readcounts.tsv**: Number of reads in each sample's FastQ after Kneaddata cleaning. 
    - **flagstat_sum.txt**: Proportion of reads in each sample's FastQ that are either uniquely mapped or chimeric over all mapped reads **to the assembled scaffolds** (the other category is multi-mapped).
    - **sr_catalogue_mapreads.tab**: Proportion of reads in each sample's FastQ that are either uniquely mapped or chimeric over all mapped reads **to the refined MAG bins** (the other category is multi-mapped).
    - **coas_catalogue_mapreads.tab**: sr_catalogue_mapreads.tab: Proportion of reads in the FastQs in each coassembly set that are either uniquely mapped or chimeric over all mapped reads **to the refined MAG bins** (the other category is multi-mapped).
    - **sr_checkm_metrics.csv**: Completeness, contamination, and strain heterogeneity of each single sample-assembled refined MAG bin.
    - **coas_checkm_metrics.csv**: Completeness, contamination, and strain heterogeneity of each coassembled refinedMAG bin.
    - **dnadiff_summary.tsv**: Comparison between each refined MAG bin (single-sample and coassembly) and the most similar RefSeq reference genome. Includes total bases in reference genome, the number of aligned bases, and the average identity. 

### Generating Figures

- The R commands have to be run separately because R also has GLIBC_2.14 issues in the Snakemake environment.
```R 
setwd("/path/to/final/reports")
source("~/bin/MAG_Snakemake_wf/scripts/plotting/plot_sc.R") # plot_ss.R for single-sample-only 

checkm_df = graph_checkm()
graph_cmseq()
graph_alignment()
graph_dnadiff(checkm_df)
graph_gtdb()
```

### Environment 

```Bash
checkm data setRoot /athena/masonlab/scratch/databases/metagenomics/CheckM_2015_01_16
export GTDBTK_DATA_PATH=/athena/masonlab/scratch/databases/metagenomics/GTDBTk_R202
export PATH="/home/lam4003/bin/MAG_Snakemake_wf/scripts:/home/lam4003/bin/MaxBin-2.2.7:/home/lam4003/bin/MaxBin-2.2.7/auxiliary/FragGeneScan_1.30:/home/lam4003/bin/MaxBin-2.2.7/auxiliary/hmmer-3.1b1/src:/home/lam4003/bin/MaxBin-2.2.7/auxiliary/bowtie2-2.2.3:/home/lam4003/bin/MaxBin-2.2.7/auxiliary/idba-1.1.3/bin:$PATH" # 1,2
export LD_LIBRARY_PATH=/home/lam4003/bin/anaconda3/envs/snakemake/glibc-2.14/lib # 3
conda activate snakemake
```
1. `preprocessing.Snakefile: kneaddata_bowtie`: BBMap internal code requires `scripts/` to be PATH-accessible.
2. `binning.Snakefile: metawrap, metawrap_coas`: For MaxBin2, which was not installed through conda.
3. `framework.Snakefile: mapreads_scaffold`: May be required.

### Local Job Submission

```Bash
snakemake --keep-going --restart-times 0 --jobs 50 --use-conda --conda-prefix /home/lam4003/bin/anaconda3/envs
```
- Mainly used to stress test the pipeline on a subset of the data, hence 0 restarts

### Scheduled Job Submission

```Bash
sbatch -J job_name ~/bin/MAG_Snakemake_wf/run_msw.sh
```
- The following Snakemake flags are set:
    - `--keep-going`: Continue with the rest of the rule script upon encountering errors.
    - `--restart-times 3`: Number of times to restart failed jobs.
    - `--jobs 50`: Maximum number of Slurm jobs in progress at once.
    - `--use-conda`: Tells the pipeline to use external conda environments.
    - `--conda-prefix`: Tells the pipeline where to find the external conda environments.
    - `--cluster`: The Slurm submission command for each rule.

## Modifications

- For all rules and `~/bin/metawrap-mg-1.3.2/bin/metawrap-modules/binning.sh` that refer to a script (ex. `*.sh` or `*.pl`), the implicit call was replaced by a hard path to minimize PATH interpretation issues.

### Input and Output

- `coassembly_runs.txt`: Each row is now a list of samples associated with a coassembled set, instead of explicit paths to the post-processed FastQs for ease of human input as well as code generalizability. 
- `clusterconfig.yaml`: Standardized job names to include sample-related wildcard information, step name, and Slurm job number.
- Instead of outputting finished graphics, final reports are gathered instead for easy downloading. 

### Preprocessing

- There is one version for direct downloads of gzipped FastQs from the SRA (`preprocessing.Snakefile`), and another for locally available data (`preprocessing_ld.Snakefile`).
- `kneaddata_bowtie`: Had to copy all of BBMap's Java dependencies to `~/bin/MAG_Snakemake_wf/scripts`
- `kneaddata_bowtie`: Hard-codes the location of the human genome index at the common GTDB-Tk directory i.e.: step is not dependent on a local `data/` directory.
- `bayeshammer`: Now a rule in preprocessing as opposed to a step attached to SPAdes. Previously, Bayeshammer was i) rerun every time SPAdes was and ii) doubly run for the coassembled samples, which wasted a ton of time and storage.
- Most of the `coas.Snakefile` rules have been eliminated because they essentially do nohing. 

### De Novo Assembly

- The following apply to `spades` and `spades_coas`
- Call `est_mem_spades` and `est_mem_coas` respectively to estimate the GB required based on the size of the dataset. Linear model is loosely based on MSW estimates, some MIAB Illumina runs (102TP, coassembly), and the Almeida test runs. 
    - Number of attempts is kept in a separate text file, instead of being calculated by Snakemake in case master job has to restart. 
- SPAdes reruns over and over despite having completed prior iterations of k! Introduced a check for finished iterations based on existing assembly graphs.
- SPAdes now takes unpaired reads generated by BayesHammer from preprocessing

### MAG Binning

- `metawrap` and `metawrap_coas`: Handles the generation of zero bins appropriately.
- `concoct_binning` and `concoct_binning_coas`: CONCOCT was made into its own rule, since it required a separate environment. It requires MetaWRAP to be finished since it depends on some of its intermediate files.
    - Not a custom module- is basically the CONCOCT step of MetaWRAP copied into its own script, `concoct.sh`.

### MAG Bin Refinement

- `bin_refinement` and `refine_bins_init`: Allows for MetaWRAP and CONCOCT to generate no bins without having the steps fail. Now depends on CONCOCT finishing.
- `dRep`: Handles the case where there are one or fewer bins appropriately. 

### Plotting

- R scripts for graphics generation were removed from 5 Snakefiles because of the `GLIBC_2.14` incompatibilities for a majority of them.
- Plotting scripts will be unified to remove redundancies, extended, and automated in the near future. 

### Quality Checking

- `dnadiff`: Doesn't take gunzipped FastAs, which are the storage format of RefSeq sequences. RefSeq FastAs are temporarily extracted for `dnadiff` and deleted upon completion.
- CheckM and GTDB do not set the path to their respective databases within the respective rules. 

## More Information

- No citations available yet for this version of the pipeline
- For a walk-through of the original pipeline and sample figures, please read and cite: Saheb Kashaf, S., Almeida, A., Segre, J.A. et al. Recovering prokaryotic genomes from host-associated, short-read shotgun metagenomic sequencing data. Nat Protoc (2021). https://doi.org/10.1038/s41596-021-00508-2
