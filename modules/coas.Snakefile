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


def get_sample_reads_coas(sample):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    print(df, file=sys.stderr)
    datasets = df[df["coassembly"] == sample]["datasets"][0].split(",") # S1,S2
    prefix = join(DATA_DIR, preprocessing_dir, "processed", "singlerun")
    r1 = [join(prefix, i + "_1.fastq.gz") for i in datasets]
    r2 = [join(prefix, i + "_2.fastq.gz") for i in datasets]
    co = [join(prefix, i + "_upd.fastq.gz") for i in datasets]
    dict = {"r1": r1, "r2": r2, "co": co}
    print(dict)
    return dict


rule coas:
    input:
        unpack(lambda wildcards: get_sample_reads_coas(wildcards.sample)),
    output:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{sample}_1.fastq.gz"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{sample}_2.fastq.gz"),
        upd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{sample}_upd.fastq.gz"),
    shell:
        """
        cat {input.r1} > {output.fwd}
        cat {input.r2} > {output.rev}
        cat {input.co} > {output.upd}
        """


# rule sort:
#     input:
#         fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_1.fastq"),
#         rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_2.fastq"),
#     output:
#         fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_1.fastq"),
#         rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_2.fastq"),
#     singularity:
#         "shub://sskashaf/MAG_wf_containers_2021:assembly"
#     shell:
#         """
#         /home/lam4003/bin/bbmap/repair.sh in={input.fwd} in2={input.rev} out={output.fwd} out2={output.rev} repair
#         rm -f {input.fwd}
#         rm -f {input.rev}
#         """


# rule zip_file:
#     input:
#         fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_1.fastq"),
#         rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_2.fastq"),
#     output:
#         fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_1.fastq.gz"),
#         rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_2.fastq.gz"),
#     params:
#         fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_1.fastq.gz"),
#         rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_2.fastq.gz"),
#     shell:
#         """
#         rm -f {params.fwd}
#         rm -f {params.rev}
#         rm -f {output.fwd}
#         rm -f {output.rev}
#         gzip {input.fwd}
#         gzip {input.rev}
#         mv {params.fwd} {output.fwd}
#         mv {params.rev} {output.rev} 
#         """




# rule readcounts_coas:
#     input:
#         raw=expand(join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_{read}.fastq"), run=COAS, read=["1", "2"]),
#     output:
#         join(DATA_DIR, preprocessing_dir, "readcounts_coas.tsv"),
#     run:
#         outfile = str(output)
#         if os.path.exists(outfile):
#             os.remove(outfile)
#         with open(outfile, "w") as outf:
#             outf.writelines("Run\tReadcount\n")
#             for run in input:
#                 readcount = int(linecount(run))
#                 line = run + "\t" + str(readcount)
#                 outf.writelines(line + "\n")
