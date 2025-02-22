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


import shutil


# READS MAPPING TO ASSEMBLY

rule mapreads_scaffold:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_2.fastq"),
        scaffold=join(DATA_DIR, assembly_dir, "singlerun/{run}/scaffolds.fasta"),
    output:
        flagstat=join(DATA_DIR, assembly_dir, "singlerun/{run}/mapreads/flagstat.txt"),
    params:
        dir=join(DATA_DIR, assembly_dir, "singlerun/{run}/mapreads"),
        scaffold=join(DATA_DIR, assembly_dir, "singlerun/{run}/{run}_scaffolds.fasta"),
        alignedsorted=join(DATA_DIR, assembly_dir, "singlerun//{run}/mapreads/alignedsorted.bam"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    shell:
        """
        rm -rf {params.dir}
        mkdir -p {params.dir}
        scp {input.scaffold} {params.scaffold}
        bwa index {params.scaffold}
        bwa mem -t {threads} {params.scaffold} {input.fwd} {input.rev} | samtools view -bS - | \
        samtools sort -@ {threads} -o {params.alignedsorted} -
        samtools index {params.alignedsorted}
        samtools flagstat {params.alignedsorted} > {output.flagstat}
        rm {params.alignedsorted}
        """


rule parse_mapreads_scaffold:
    input:
        join(DATA_DIR, assembly_dir, "singlerun/{run}/mapreads/flagstat.txt"),
    output:
        join(DATA_DIR, assembly_dir, "singlerun/{run}/mapreads/flagstat_parsed.txt"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    shell:
        """
        crimson flagstat {input}>{output}
        """


rule aggregate_mapreads_scaffold:
    input:
        expand(join(DATA_DIR, assembly_dir, "singlerun/{run}/mapreads/flagstat_parsed.txt"), run=RUN),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/flagstat_sum.txt"),
    run:
        with open(str(output), "w") as outf:
            for i in input:
                run = str(i).split("singlerun/")[1]
                run = run.split("/mapreads/")[0]
                dict = eval(open(str(i), "r").read())
                supplementary = dict["pass_qc"]["supplementary"]
                secondary = dict["pass_qc"]["secondary"]
                mapped = dict["pass_qc"]["mapped"]
                total = dict["pass_qc"]["total"]
                perc = (mapped - supplementary - secondary) 
                outf.write(str(run) + "\t" + str(perc) + "\n")
        outf.close()





# READ MAPPING TO MAGS

rule cat_MAGs:
    input:
        out=join(DATA_DIR, binning_analyses, "singlerun/dRep/done.txt"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/framework/bwa-ref_name_vf/ref-db.fasta"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    params:
        indir=join(DATA_DIR, binning_analyses, "singlerun/dRep/dereplicated_genomes/"),
    shell:
        """
        cat {params.indir}/*fa>{output}
        bwa index {output}
        """


rule cat_MAGs_coas:
    input:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/dRep/data_tables/Sdb.csv"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/bwa-ref_name_vf/ref-db.fasta"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    params:
        indir=join(DATA_DIR, binning_analyses, "singlerun_coassembly/dRep/dereplicated_genomes/"),
    shell:
        """
        cat {params.indir}/*.fa>{output}
        bwa index {output}
        """


rule readmap:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_2.fastq"),
        catalogue=join(DATA_DIR, binning_analyses, "singlerun/framework/bwa-ref_name_vf/ref-db.fasta"),
    output:
        flagstat=join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/flagstat/{run}.txt"),
    params:
        alignedsorted=join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/flagstat/tmp_{run}.bam"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    shell:
        """
        bwa mem -t {threads} {input.catalogue} {input.fwd} {input.rev} \
        | samtools view -bS - | samtools sort -@ {threads} -o {params.alignedsorted} -
        samtools index {params.alignedsorted}
        samtools flagstat {params.alignedsorted} > {output.flagstat}
        rm {params.alignedsorted}
        """


rule readmap_coassembly:
    input:
        # fwd=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_1.fastq"),
        # rev=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun/{run}_2.fastq"),
        catalogue=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/bwa-ref_name_vf/ref-db.fasta"),
    output:
        flagstat=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/flagstat/{run}.txt"),
    params:
        kdb_dir=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/singlerun"),
        alignedsorted=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/flagstat/tmp_{run}.bam"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    shell:
        """
        cat {params.kdb_dir}/*_1.fastq > {params.kdb_dir}/tmp_1.fastq
        cat {params.kdb_dir}/*_2.fastq > {params.kdb_dir}/tmp_2.fastq
        bwa mem -t {threads} {input.catalogue} {params.kdb_dir}/tmp_1.fastq {params.kdb_dir}/tmp_2.fastq \
        | samtools view -bS - | samtools sort -@ {threads} -o {params.alignedsorted} -
        samtools index {params.alignedsorted}
        samtools flagstat {params.alignedsorted} > {output.flagstat}
        rm {params.alignedsorted}
        rm {params.kdb_dir}/tmp*
        """


rule parse_readmap:
    input:
        join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/flagstat/{run}.txt"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/flagstat_parsed/{run}.txt"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    shell:
        """
        crimson flagstat {input}>{output}
        """


rule parse_readmap_coas:
    input:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/flagstat/{run}.txt"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/flagstat_parsed/{run}.txt"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    shell:
        """
        crimson flagstat {input}>{output}
        """


rule write_scaffold:
    input:
        expand(join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/flagstat_parsed/{run}.txt"), run=RUN),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/sr_catalogue_mapreads.tab"),
    run:
        with open(str(output), "w") as outf:
            for i in input:
                base = os.path.basename(str(i))
                run = os.path.splitext(base)[0]
                dict = eval(open(str(i), "r").read())
                supplementary = dict["pass_qc"]["supplementary"]
                secondary = dict["pass_qc"]["secondary"]
                mapped = dict["pass_qc"]["mapped"]
                total = dict["pass_qc"]["total"]
                perc = (mapped - supplementary - secondary) 
                outf.write(str(run) + "\t" + str(perc) + "\n")
        outf.close()



rule write_scaffold_coas:
    input:
        expand(join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/flagstat_parsed/{run}.txt"), run=COAS),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/coas_catalogue_mapreads.tab"),
    run:
        with open(str(output), "w") as outf:
            for i in input:
                base = os.path.basename(str(i))
                run = os.path.splitext(base)[0]
                dict = eval(open(str(i), "r").read())
                supplementary = dict["pass_qc"]["supplementary"]
                secondary = dict["pass_qc"]["secondary"]
                mapped = dict["pass_qc"]["mapped"]
                total = dict["pass_qc"]["total"]
                perc = (mapped - supplementary - secondary) 
                outf.write(str(run) + "\t" + str(perc) + "\n")
        outf.close()


# rule plot_framework:
#     input:
#         mapreads_sr=join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/sr_catalogue_mapreads.tab"),
#         mapreads_coas=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/coas_catalogue_mapreads.tab"),
#         flagstat=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/flagstat_sum.txt"),
#         readcounts=join(DATA_DIR, preprocessing_dir, "readcounts.tsv"),
#         dRep=join(DATA_DIR, binning_analyses, "singlerun/dRep/done.txt"),
#     output:
#         join(DATA_DIR, "figures/perassemb_perref.png"),
#     singularity:
#         "shub://sskashaf/MAG_wf_containers_2021:r"
#     params:
#         summary=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/summary_framework.csv")
#     shell:
#         """
#         Rscript /home/lam4003/bin/MAG_Snakemake_wf/scripts/plotting/plot_framework.R {input.readcounts} {input.flagstat} {input.mapreads_sr} {input.mapreads_coas} {params.summary}
#         """




# ASSEMBLY AND BINNING SUMMARY STATISTICS


rule make_reference_dir:
    input:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/GTDB/gtdbtk.bac120.summary.tsv"),
    output:
        join(DATA_DIR, assembly_dir, "metaQUAST/done_sc.txt"),
    params:
        refdir=join(DATA_DIR, assembly_dir, "metaQUAST/references"),
    run:
        df = pd.read_csv(str(input), sep = '\t')
        refs_tmp = list(df.iloc[:,2].unique())
        refs = [r for r in refs_tmp if str(r) != 'nan']
        GTDB = os.getenv('GTDBTK_DATA_PATH')
        if not os.path.isdir(params.refdir): os.makedirs(params.refdir)
        for r in refs: # GCF_000190535.1
            parts = str(r).split('_')
            path = join(GTDB, 'fastani/database', parts[0], parts[1][0:3], parts[1][3:6], parts[1][6:9], r + '_genomic.fna.gz')
            shutil.copy(path, params.refdir)
        open(str(output), "w").close()


rule make_reference_dir_ss:
    input:
        join(DATA_DIR, binning_analyses, "singlerun/GTDB/gtdbtk.bac120.summary.tsv"),
    output:
        join(DATA_DIR, assembly_dir, "metaQUAST/done_ss.txt"),
    params:
        refdir=join(DATA_DIR, assembly_dir, "metaQUAST/references"),
    run:
        df = pd.read_csv(str(input), sep = '\t')
        refs_tmp = list(df.iloc[:,2].unique())
        refs = [r for r in refs_tmp if str(r) != 'nan']
        GTDB = os.getenv('GTDBTK_DATA_PATH')
        if not os.path.isdir(params.refdir): os.makedirs(params.refdir)
        for r in refs: # GCF_000190535.1
            parts = r.split('_')
            path = join(GTDB, 'fastani/database', parts[0], parts[1][0:3], parts[1][3:6], parts[1][6:9], r + '_genomic.fna.gz')
            shutil.copy(path, params.refdir)
        open(str(output), "w").close()


def find_assembly_strategy(wildcards):
    if open("Snakefile", "r").readline().strip().split()[1] == "COASSEMBLY":
        return "data/01_assembly/metaQUAST/done_sc.txt"
    else:
        return "data/01_assembly/metaQUAST/done_ss.txt"

rule metaquast:
    input:
        lambda wildcards: find_assembly_strategy(wildcards), 
    output:
        join(DATA_DIR, assembly_dir, "metaQUAST/combined_reference/report.tsv"),
    params:
        outdir=join(DATA_DIR, assembly_dir, "metaQUAST"),
        refdir=join(DATA_DIR, assembly_dir, "metaQUAST/references"),
        mincontiglength=1000, # Report on only the contigs that are binned
    threads: workflow.cores,
    shell:
        """
        scaffolds=`ls data/01_assembly/*/*/scaffolds.fasta`
        python /home/lam4003/bin/quast-master/metaquast.py --threads {threads} -r {params.refdir} -m {params.mincontiglength} -o {params.outdir} $scaffolds --no-plots
        """

# rule write_binning_stats:
#     input:
#         join(DATA_DIR, binning_dir, "singlerun_binning.csv"),
#         join(DATA_DIR, binning_dir, "coassembly_binning.csv"),
#     output:
#         join(DATA_DIR, binning_dir, "binning_stats.csv"), # Not combined_reference/ because no reference genomes
#     shell:
#         """
#         cat {input} > {output}
#         """


