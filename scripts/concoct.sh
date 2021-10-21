########################################################################################################
########################                    RUNNING CONCOCT                     ########################
########################################################################################################

error () { echo "$1" "*"; exit 1; }

echo "RUNNING CONCOCT"

if [[ ! -s ${1}/work_files/concoct_depth.txt ]]; then 
    echo "indexing .bam alignment files..."
    for FILE in ${1}/work_files/*.bam; do
        echo $FILE
        samtools index -@ $threads -b $FILE
    done

    echo "cutting up contigs into 10kb fragments for CONCOCT..."
    /home/lam4003/bin/anaconda3/envs/binning/bin/cut_up_fasta.py ${1}/work_files/assembly.fa -c 10000 --merge_last -b ${1}/work_files/assembly_10K.bed -o 0 > ${1}/work_files/assembly_10K.fa 
        if [[ $? -ne 0 ]]; then error "Something went wrong with cutting up contigs. Exiting."; fi

    echo "estimating contig fragment coverage..."   
    CMD="/home/lam4003/bin/anaconda3/envs/binning/bin/concoct_coverage_table.py ${1}/work_files/assembly_10K.bed ${1}/work_files/*.bam > ${1}/work_files/concoct_depth.txt" 
    $(eval $CMD)
    if [[ $? -ne 0 ]]; then error "Something went wrong with estimating fragment abundance. Exiting..."; fi
else
    echo "looks like contig coverage was already estimated... skipping"
fi

echo "Starting binning with CONCOCT..."
mkdir ${1}/work_files/concoct_out

concoct -l $2 \
    --coverage_file ${1}/work_files/concoct_depth.txt \
    --composition_file ${1}/work_files/assembly_10K.fa \
    -b ${1}/work_files/concoct_out

if [[ $? -ne 0 ]]; then error "Something went wrong with binning with CONCOCT. Exiting..."; fi
sed  -i '1i contig_id,cluster_id' ${1}/work_files/concoct_out/clustering_gt${2}.csv

echo "merging 10kb fragments back into contigs"
/home/lam4003/bin/anaconda3/envs/binning/bin/merge_cutup_clustering.py ${1}/work_files/concoct_out/clustering_gt${2}.csv > ${1}/work_files/concoct_out/clustering_gt${2}_merged.csv  
if [[ $? -ne 0 ]]; then error "Something went wrong with merging fragments. Exiting..."; fi

echo "splitting contigs into bins"
mkdir ${1}/concoct_bins
/home/lam4003/bin/metawrap-mg-1.3.2/bin/metawrap-scripts/split_concoct_bins.py ${1}/work_files/concoct_out/clustering_gt${2}_merged.csv ${1}/work_files/assembly.fa ${1}/concoct_bins 
echo "CONCOCT finished successfully, and found $(ls -l ${1}/concoct_bins | grep .fa | wc -l) bins!"