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


def get_sample_reads_cmseq(sample):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    print(df, file=sys.stderr)
    datasets = df[df["coassembly"] == sample]["datasets"][0].split(",") # S1,S2
    prefix = join(DATA_DIR, preprocessing_dir, "singlerun")
    return [join(prefix, i + ".fastq.gz") for i in datasets]


rule rename_fasta:
    input:
        MAG=join(DATA_DIR, binning_dir, "singlerun/{sample}/refined_bins_50_10/{file}.fa"),
    output:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
    params:
        name="{sample}_{file}",
    shell:
        """
        /home/lam4003/bin/MAG_Snakemake_wf/scripts/rename_multifasta_prefix.py -f {input.MAG} -p {params.name} > {output}
        """


rule prokka:
    input:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
    output:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.gff",
        ),
    conda:
        "/home/lam4003/bin/anaconda3/envs/prokka.yaml"
    singularity:
        "docker://staphb/prokka:1.14.5"
    params:
        out_prokka=join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/"),
        prefix="{sample}_metawrap_refined_{file}",
    shell:
        """
        prokka {input} --kingdom Bacteria --outdir {params.out_prokka} \
        --prefix {params.prefix} --force --locustag {params.prefix} --cpus {threads}
        """


rule rename_fasta_coas:
    input:
        MAG=join(DATA_DIR, binning_dir, "coassembly/{sample}/refined_bins_50_10/{file}.fa"),
    output:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
    params:
        name="{sample}_{file}",
    shell:
        """
        /home/lam4003/bin/MAG_Snakemake_wf/scripts/rename_multifasta_prefix.py -f {input.MAG} -p {params.name} > {output}
        """


rule prokka_coas:
    input:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
    output:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.gff",
        ),
    conda:
        "/home/lam4003/bin/anaconda3/envs/prokka.yaml"
    singularity:
        "docker://staphb/prokka:1.14.5"
    params:
        out_prokka=join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/"),
        prefix="{sample}_metawrap_refined_{file}",
    shell:
        """
        prokka {input} --kingdom Bacteria --outdir {params.out_prokka} \
        --prefix {params.prefix} --force --locustag {params.prefix} --cpus {threads}
        """


rule cmseq:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{sample}.fastq.gz"),
        MAG=join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
        prokka=join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.gff",
        ),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.cmseq.csv"),
    threads: workflow.cores
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:cmseq"
    params:
        name=join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}"),
        dir=join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/"),
    shell:
        """
        /home/lam4003/bin/MAG_Snakemake_wf/scripts/cmseq.sh -t {threads} -i {input.fq} -r {input.MAG} -g {input.prokka} -o {params.name}
        """


rule cmseq_coas:
    input:
        MAG=join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
        prokka=join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.gff",
        ),
    output:
        cmseq=join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.cmseq.csv",
        ),
        done=join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_{file}_done.txt",
        ),
    threads: workflow.cores
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:cmseq"
    params:
        fq=lambda wildcards: get_sample_reads_cmseq(wildcards.sample),
        name=join(
            DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}"
        ),
    shell:
        """
        rm -f {output.cmseq}
        rm -f {output.done}
        for i in {params.fq}; do run=$(basename ${{i}} .fastq.gz); /home/lam4003/bin/MAG_Snakemake_wf/scripts/cmseq.sh\
        -t {threads} -i ${{i}} -r {input.MAG} -g {input.prokka}\
        -o {params.name}_${{run}}; awk 'NR==2' {params.name}_${{run}}.cmseq.csv >> {output.cmseq};done

        touch {output.done}
        """


rule append_cmseq:
    input:
        lambda wildcards: expand(rules.cmseq.output, file=aggregate_bins(wildcards), sample=wildcards.sample),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/sum_cmseq_{sample}.csv"),
    params:
        join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/"),
    run:
        if os.path.exists(str(output)):
            os.remove(str(output))
        else:
            print("Aggregating single run strain heterogeneity results")
        path = str(params)
        file1 = open(str(output), "w")
        for filename in glob.glob(os.path.join(path, "*.csv")):
            with open(filename, "r") as f:
                lines = f.readlines()
                if len(lines) > 1:
                    lines_sub = lines[1].strip()
                    L = filename + "\t" + lines_sub + "\n"
                    file1.writelines(L)
        file1.close()



rule append_cmseq_coas:
    input:
        lambda wildcards: expand(rules.cmseq_coas.output, file=aggregate_bins_coas(wildcards), sample=wildcards.sample),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/sum_cmseq_{sample}.csv"),
    params:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/"),
    run:
        if os.path.exists(str(output)):
            os.remove(str(output))
        else:
            print("Aggregating coassembly strain heterogeneity results")
        path = str(params)
        file1 = open(str(output), "w")
        for filename in glob.glob(os.path.join(path, "*SRR*.csv")):
            with open(filename, "r") as f:
                lines = f.readlines()
                if len(lines) > 1:
                    lines_sub = lines[1].strip()
                    L = filename + "\t" + lines_sub + "\n"
                    file1.writelines(L)
        file1.close()



rule aggregate_cmseq:
    input:
        sr=expand(join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/sum_cmseq_{sample}.csv"), sample=RUN),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/cmseq/summ_cmseq_all.csv"),
    shell:
        """
        cat {input}>{output}         
        """


rule aggregate_cmseq_coas:
    input:
        sr=expand(join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/sum_cmseq_{sample}.csv"), sample=COAS),
    output:
        coas=join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/summ_cmseq_all.csv"),
    shell:
        """
        cat {input}>{output}
        """


# rule plot_cmseq:
#     input:
#         coas=join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/summ_cmseq_all.csv"),
#         sr=join(DATA_DIR, binning_analyses, "singlerun/cmseq/summ_cmseq_all.csv"),
#     output:
#         join(DATA_DIR, "figures/cmseq_plot.png"),
#     singularity:
#         "shub://sskashaf/MAG_wf_containers_2021:r"
#     shell:
#         """
# 	if [ -s {input.sr} ] && [ -s {input.coas} ];
#         then
#                 Rscript /home/lam4003/bin/MAG_Snakemake_wf/scripts/plotting/plot_cmseq.R {input.sr} {input.coas}
#         else
#                 echo "{input.sr} or {input.coas} are empty. Empty plot printed."
#                 touch {output} # Creates an empty plot file to satisfy requirements
#         fi
#         """
