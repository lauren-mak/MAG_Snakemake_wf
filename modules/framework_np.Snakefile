
# READS MAPPING TO ASSEMBLY

# Replaced bwa (short reads) with minimap2 (long reads)
rule mapreads_scaffold:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
        asm=join(DATA_DIR, assembly_dir, "final_polished/singlerun/{run}.polished.fasta"),
    output:
        flagstat=join(DATA_DIR, assembly_dir, "singlerun/{run}/mapreads/flagstat.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    params:
        dir=join(DATA_DIR, assembly_dir, "singlerun/{run}/mapreads"),
        asm=join(DATA_DIR, assembly_dir, "singlerun/{run}/{run}_scaffolds.fasta"),
        alignedsorted=join(DATA_DIR, assembly_dir, "singlerun//{run}/mapreads/alignedsorted.bam"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    shell:
        """
        rm -rf {params.dir}
        mkdir -p {params.dir}
        scp {input.asm} {params.asm}
        minimap2 -ax map-ont {input.asm} {input.fq} | samtools view -bS - | samtools sort -@ {threads} -o {params.alignedsorted} -
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
        """


rule readmap:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "singlerun/{run}.fastq.gz"),
        catalogue=join(DATA_DIR, binning_analyses, "singlerun/framework/bwa-ref_name_vf/ref-db.fasta"),
    output:
        flagstat=join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/flagstat/{run}.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    params:
        alignedsorted=join(DATA_DIR, binning_analyses, "singlerun/framework/mapreads/flagstat/tmp_{run}.bam"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    shell:
        """
        minimap2 -ax map-ont {input.catalogue} {input.fq} | samtools view -bS - | samtools sort -@ {threads} -o {params.alignedsorted} -
        samtools index {params.alignedsorted}
        samtools flagstat {params.alignedsorted} > {output.flagstat}
        rm {params.alignedsorted}
        """


rule readmap_coassembly:
    input:
        fq=join(DATA_DIR, preprocessing_dir, "coassembly/{run}.fastq.gz"),
        catalogue=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/bwa-ref_name_vf/ref-db.fasta"),
    output:
        flagstat=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/flagstat/{run}.txt"),
    conda:
        "/home/lam4003/bin/anaconda3/envs/nanopore.yaml"
    params:
        alignedsorted=join(DATA_DIR, binning_analyses, "singlerun_coassembly/framework/mapreads/flagstat/tmp_{run}.bam"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:framework"
    shell:
        """
        minimap2 -ax map-ont {input.catalogue} {input.fq} | samtools view -bS - | samtools sort -@ {threads} -o {params.alignedsorted} -
        samtools index {params.alignedsorted}
        samtools flagstat {params.alignedsorted} > {output.flagstat}
        rm {params.alignedsorted}
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
        refs = list(df.iloc[:,2].unique())
        if "N/A" in refs: refs = refs.remove("N/A")
        GTDB = os.getenv('GTDBTK_DATA_PATH')
        os.mkdirs(params.refdir)
        for r in refs: # GCF_000190535.1
            parts = r.split('_')
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
        refs = list(df.iloc[:,2].unique())
        if "N/A" in refs: refs = refs.remove("N/A")
        GTDB = os.getenv('GTDBTK_DATA_PATH')
        os.mkdirs(params.refdir)
        for r in refs: # GCF_000190535.1
            parts = r.split('_')
            path = join(GTDB, 'fastani/database', parts[0], parts[1][0:3], parts[1][3:6], parts[1][6:9], r + '_genomic.fna.gz')
            shutil.copy(path, params.refdir)
        open(str(output), "w").close()


def find_assembly_strategy(sample):
    if open("Snakefile", "r").readline().strip().split()[1] == "COASSEMBLY":
        return "data/01_assembly/metaQUAST/done_sc.txt"
    else:
        return "data/01_assembly/metaQUAST/done_ss.txt"

rule metaquast:
    input:
        lambda wildcards: find_assembly_strategy(wildcards), 
    output:
        join(DATA_DIR, assembly_dir, "metaQUAST/report.tsv"), # Not combined_reference/ because no reference genomes
    params:
        outdir=join(DATA_DIR, assembly_dir, "metaQUAST"),
        mincontiglength=1000, # Report on only the contigs that are binned
    shell:
        """
        scaffolds=`ls data/01_assembly/final_polished/*/*.polished.fasta`
        python /home/lam4003/bin/quast-master/metaquast.py --threads {threads} --max-ref-number 0 -m {params.mincontiglength} -o {params.outdir} $scaffolds
        """

