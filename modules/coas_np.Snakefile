# def get_nano_reads(sample):
#     sample_reads = []
#     sample_file = "coassembly_runs.txt"
#     df = pd.read_csv(sample_file, sep="\t")
#     print(df, file=sys.stderr)
#     datasets = df[df["coassembly"] == sample]["datasets"][0].split(",") # S1,S2
#     prefix = join(DATA_DIR, preprocessing_dir, "filtlong")
#     return [join(prefix, i + ".fastq.gz") for i in datasets]

# rule merge:
#     input:
#         unpack(lambda wildcards: get_nano_reads(wildcards.sample)),
#     output:
#         join(DATA_DIR, preprocessing_dir, "processed/{run}.fastq.gz"),
#     params:
#         short_qc="2500",
#     threads: workflow.cores
#     shell:
#         """
#         cat {input} > {output}
#         """
