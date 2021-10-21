#!/bin/bash
#SBATCH --partition=panda

snakemake --keep-going --restart-times 3 --jobs 50 --use-conda --conda-prefix /home/lam4003/bin/anaconda3/envs --cluster-config clusterconfig.yaml --cluster "sbatch -n {cluster.nCPU} --mem {cluster.mem} -e {cluster.error} -o {cluster.output} -t {cluster.time}"

