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

# concoct=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/concoct_bins"),
# -C {params.concoct}
# checkm data setRoot ~/checkm_database/
rule bin_refinement_coas:
    input:
        join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/concoct_done.txt"),
    output:
        join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/metawrap_50_10_bins/done.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/binning.yaml"
    singularity:
        "shub://sskashaf/Containers:metawrap"
    params:
        metabat=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/metabat2_bins"),
        maxbin=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/maxbin2_bins"),  
        concoct=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/concoct_bins"),
        outdir=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement"),
    threads: workflow.cores
    shell:
        """ 
        rm -rf {params.outdir}
        metawrap bin_refinement -o {params.outdir} -t {threads} -A {params.metabat} -B {params.maxbin} -C {params.concoct} -c 50 -x 10
        touch {output}
        """


checkpoint refine_bins_init_coas:
    input:
        rules.bin_refinement_coas.output,
    output:
        directory(join(DATA_DIR, binning_dir, "coassembly/{sample}/refined_bins_50_10")),
    params:
        refined_folder=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/metawrap_50_10_bins"),
        folder=join(DATA_DIR, binning_dir, "coassembly/{sample}/refined_bins_50_10/"),
    shell:
        """
        rm -rf {params.folder}
        mkdir -p {params.folder}
        cp {params.refined_folder}/*.fa {params.folder}
        """


def aggregate_bins_coas(wildcards):
    refine_out = checkpoints.refine_bins_init_coas.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(refine_out, "{bin}.fa")).bin


rule copy_bins_coas:
    input:
        join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/metawrap_50_10_bins/{i}.fa"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/all_metawrap_bins_coas/{sample}_metawrap_refined_{i}.fa"),
        join(DATA_DIR, "all_bins/{sample}_metawrap_refined_{i}.fa"),
    params:
        dir1=join(DATA_DIR, binning_analyses, "singlerun_coassembly/all_metawrap_bins_coas/"),
        dir2=join(DATA_DIR, "all_bins/"),
        out1=join(DATA_DIR, binning_analyses, "singlerun_coassembly/all_metawrap_bins_coas/{sample}_metawrap_refined_{i}.fa"),
        out2=join(DATA_DIR, "all_bins/{sample}_metawrap_refined_{i}.fa"),
    shell:
        """
        mkdir -p {params.dir1}
        mkdir -p {params.dir2}
        cp {input} {params.out1}
        cp {input} {params.out2}
        """


rule copy_refined_coas:
    input:
        lambda wildcards: expand(rules.copy_bins_coas.output, i=aggregate_bins_coas(wildcards), sample=wildcards.sample),
    output:
        join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/copied.txt"),
    shell:
        """
        touch {output}
        """

# checkm data setRoot ~/checkm_database/
rule checkm_coas:
    input:
        expand(join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/copied.txt"), sample=COAS),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.tsv"),
    params:
        ext="fa",
        indir=join(DATA_DIR, binning_analyses, "singlerun_coassembly/all_metawrap_bins_coas"),
        outdir=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/"),
    singularity:
        "shub://sskashaf/Containers:metawrap"
    shell:
        """ 
        rm -rf {params.outdir}
        checkm lineage_wf -t {threads} -x {params.ext} --tab_table -f {output} {params.indir} {params.outdir}
        """


rule parse_checkm_coas:
    input:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.tsv"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.csv"),
    params:
        checkm=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics_tmp.csv"),
        checkm2=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics_tmp2.csv"),
    shell:
        """
        sed -i '1d' {input}
        cut -f1,12,13,14 {input} | tr '\t' ','>{params.checkm}
        sed 's/,/.fa,/' {params.checkm}>{params.checkm2}
        echo -e "genome,completeness,contamination,strain_heterogeneity" | cat - {params.checkm2} > {output}
        rm {params.checkm}
        rm {params.checkm2}
        """

# rule binning_stats_coas:
#     input:
#         join(DATA_DIR, binning_dir, "singlerun_coassembly/{sample}/metawrap_bin_refinement/metawrap_50_10_bins.stats"),        
#     output: 
#         join(DATA_DIR, binning_dir, "singlerun_coassembly/{sample}.csv"),
#     params:
#         sample = "{sample}",
#     run:
#         df = pd.read_table(str(input), header = 0, sep = '\t')
#         df['Sample'] = [params['sample'] for i in range(df.shape[0])]
#         df.to_csv(output[0], index = False, sep = ',')


# def aggregate_stats_coas(wildcards):
#     s = glob_wildcards(join(DATA_DIR, binning_dir, "singlerun_coassembly/{sample}.csv"))
#     return expand(join(DATA_DIR, binning_dir, "singlerun_coassembly/{sample}.csv"), sample=s)


# rule concat_stats_coas:
#     input:
#         aggregate_stats_coas,
#     output: 
#         join(DATA_DIR, binning_dir, "coassembly_binning.csv")
#     run:
#         """
#         cat {input} > {output}
#         """


