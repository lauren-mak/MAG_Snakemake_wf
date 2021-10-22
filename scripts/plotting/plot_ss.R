library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyr)

# Figure 3a,b
# summ_cmseq_all.csv: CMSeq (strain heterogeneity) summary statistics of single-sample MAGs. 
# # Figure 3c
# checkm_metrics.csv: Completeness, contamination, and strain heterogeneity of each single sample-assembled refined MAG bin.
# Figure 3d,e
# dnadiff_summary.tsv: Comparison between each refined MAG bin (single-sample) and the most similar RefSeq reference genome. Includes total bases in reference genome, the number of aligned bases, and the average identity. 
# Figure 3f
# gtdbtk.bac120.summary.tsv: GTDB classifications of each classifiable refined MAG bin, with the ANI to the most similar RefSeq reference genome. 
# Figure 4a
# readcounts.tsv: Number of reads in each sample's FastQ after Kneaddata cleaning. 
# flagstat_sum.txt: Proportion of reads in each sample's FastQ that are either uniquely mapped or chimeric over all mapped reads to the assembled scaffolds (the other category is multi-mapped).
# catalogue_mapreads.tab: Proportion of reads in each sample's FastQ that are either uniquely mapped or chimeric over all mapped reads to the refined MAG bins (the other category is multi-mapped).

# Figure 3a,b: Box-and-whisker plot comparing completeness and contamination.
graph_checkm <- function(){
    # Load dataframes and format column names/variable types
    checkm_sr = read.csv("checkm_metrics.csv", header = TRUE, stringsAsFactor = FALSE)
    checkm_sr$status="Single-sample"

    # Make figures
    ggplot(checkm, aes(x = status, y = completeness, fill = status)) +
        geom_boxplot(outlier.shape=NA) + theme_classic() + 
        ylab("Completeness") + xlab("Assembly Approach") +
        scale_fill_manual(breaks=c("Single-sample","Coassembly"), values=c("#BBBBBB","#4477AA"))+ 
        guides(color=guide_legend(title="Approach")) + theme(legend.position="none")
    ggsave("checkm_completeness.png", width = 10, height = 10)

    ggplot(checkm, aes(x = status, y = contamination, fill = status)) +
        geom_boxplot(outlier.shape=NA) + theme_classic() +
        ylab("Contamination") + xlab("Assembly Approach") +
        scale_fill_manual(breaks=c("Single-sample","Coassembly"), values=c("#BBBBBB","#4477AA")) +
        guides(color=guide_legend(title="Approach")) + theme(legend.position="none") +
        ylim(0,10)
    ggsave("checkm_contamination.png", width = 10, height = 10)

    # Annotating MAG quality for Figure 3B (dnadiff)
    checkm$Quality = "Medium quality"
    checkm$Quality[checkm$completeness >= 90 & checkm$contamination <= 5] = "High quality"
    checkm$genome = gsub(".fa","",checkm$genome)
    return(checkm) # Required for the DNAdiff figure
} 

# Figure 3c: Box-and-whisker plot comparing strain heterogeneity.
graph_cmseq <- function() {
    # Load dataframes and format column names/variable types
    cmseq = read.csv("summ_cmseq_all.csv", header = TRUE, stringsAsFactor = FALSE)
    cmseq$status="Single-sample"
    colnames(cmseq)=c("path", "strain_het")
    cmseq$strain_het=as.numeric(as.vector(cmseq$strain_het))
    cmseq=cmseq[complete.cases(cmseq$strain_het),]

    # Make figure
    ggplot(cmseq, aes(x = status, y = strain_het, fill = status)) +
        geom_boxplot(outlier.shape=NA) + theme_classic() +
        ylab("Strain Heterogeneity") + xlab("Assembly Approach") +
        scale_fill_manual(breaks=c("Single run","Co-assembly"), values=c("#BBBBBB","#4477AA")) +
        guides(color=guide_legend(title="Approach")) + theme(legend.position="none") +
        ylim(0,12)
    ggsave("cmseq.png", width = 10, height = 10)
}

# Figure 3e: Histogram of average nucleotide identity between MAGs and identified reference genomes.
graph_dnadiff <- function(checkm_df) {
    # Load dataframes and format column names/variable types
    dnadiff = read.table("dnadiff_summary.tsv", header = TRUE, stringsAsFactor = FALSE)
    colnames(dnadiff) = c("Ref", "genome", "Ref_length", "Refcovered", "Query_length", "Queryaligned", "ANI")
    dnadiff$genome = gsub(".*/","", dnadiff$genome)
    dnadiff$genome = gsub(".fa","", dnadiff$genome)

    # Combine CheckM and DNADiff tables
    ddiff_checkm = merge(dnadiff, checkm_df, by = "genome")

    # Make figure
    ddiff_left <- ggplot(ddiff_checkm, aes(x = Queryaligned, y = Refcovered, color = status, shape = Quality)) +
        geom_point(stroke = 1,size=1) + theme_classic() +
        xlab("MAG Aligned (%)") + ylab("Reference Aligned (%)") +
        scale_color_manual(breaks=c("Single-sample", "Coassembly"), values=c("#BBBBBB", "#4477AA")) + 
        scale_shape_manual(values = c(1,2)) +
        guides(color=guide_legend(title="Approach")) +
        xlim(0,100) + ylim(0,100)
    ddiff_right <- ggplot(ddiff_checkm, aes(x = ANI, color = status)) + 
        geom_density() + theme_classic() +
        xlab("ANI") + ylab("Density") + 
        scale_color_manual(breaks=c("Single-sample","Co-assembly"), values=c("#BBBBBB", "#4477AA")) +
        guides(color=guide_legend(title="Approach")) + 
        xlim(0,100)
    ddiff_full <- grid.arrange(ddiff_left, ddiff_right, nrow = 1)
    ggsave("dnadiff.png", ddiff_full, width = 20, height = 10)
}

# Figure 3f: Stacked bar plot of single-sample and coassembled MAG taxonomic classifications.
graph_gtdb <- function() { 
    ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    # Load dataframes and filter for bacteria only
    gtdb = read.delim(paste("gtdbtk.bac120.summary.tsv", sep = ''), header = TRUE, stringsAsFactor = FALSE)
    gtdb = separate(data = gtdb, col = classification, sep = ";", into = ranks)
    bact = gtdb[which(gtdb$Domain == "d__Bacteria"),]

    # For each rank, calculate the proportions of taxon the MAGs are classified as 
    for (r in ranks[2:(length(ranks)-1)]){
        total = nrow(bact)   # total number of genomes
        rank.present = table(bact[!grepl("__$", bact[,r]),r]) # summary of genus
        total.rank = sum(rank.present) # sum of genus counts 
        if (length(rank.present) < 7){ # Take the 6 most common taxa in that rank
            start_idx = 1
        } else {
            start_idx = length(rank.present)-7
        }
        gtdb.ranks = data.frame(sort(rank.present)[start_idx:length(rank.present)])
        gtdb.ranks$Prop = gtdb.ranks$Freq/total*100
        other.name = paste(tolower(substr(r, 1, 1)), "__Other", sep="")
        novel.name = paste(tolower(substr(r, 1, 1)), "__Novel", sep="")

        other.freq = total.rank-sum(gtdb.ranks$Freq)
        novel.freq=total-total.rank  # novel genomes i.e.: unclassified MAGs

        other.prop = other.freq/total*100
        novel.prop = novel.freq/total*100

        other.ranks = data.frame(Var1=other.name, Freq=other.freq, Prop=other.prop)
        novel.ranks = data.frame(Var1=novel.name, Freq=novel.freq, Prop=novel.prop)

        ranks.fi = rbind(other.ranks, gtdb.ranks)
        ranks.fi = rbind(novel.ranks, ranks.fi)

        ranks.fi$Rank = r
        col_options = c("#999999","#CAB2D6","#A6CEE3","#FF7F00","#FB9A99","#E31A1C","#B2DF8A","#33A02C","#1F78B4")
        ranks.fi$Colour = col_options[1:nrow(ranks.fi)]

        if (exists("gtdb.fi.bac")){
            gtdb.fi.bac = rbind(gtdb.fi.bac, ranks.fi)
        } else {
            gtdb.fi.bac = ranks.fi
        }
    }

    # Make figure
    gtdb.fi.bac$Rank=factor(gtdb.fi.bac$Rank,levels=c("Phylum","Class","Order","Family","Genus"))
    tax_class <- print(ggplot(gtdb.fi.bac, aes(x=Rank, y=Prop, fill=Var1)) +
          geom_bar(stat="identity", colour="darkgrey", alpha=0.5, size=0.2, width=0.7) + theme_bw() + coord_flip() +
          ylab("Proportion of species (%)") + 
          scale_fill_manual(values=as.vector(gtdb.fi.bac$Colour), name="Taxa") + 
          scale_x_discrete(limits=ranks[(length(ranks)-1):2]) + 
          theme(legend.position="right") + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(axis.title.x = element_text(size=14)) + 
          theme(axis.text.y = element_text(size=12)) +
          theme(axis.title.y = element_blank()) +
          theme(axis.text.x = element_text(size=12)))
    ggsave("gtdb_classification.png", tax_class, width = 15, height = 10)
}

# Figure 4a: Dot plot of reads aligned to the identified reference genomes vs. binned MAGs.
graph_alignment <- function() {
    # Load dataframes and format column names/variable types
    read_counts = read.table("readcounts.tsv", header = TRUE, stringsAsFactor = FALSE)
    alned_ref = read.table("flagstat_sum.txt", header = FALSE, stringsAsFactor = FALSE)
    alned_mag_sr = read.table("catalogue_mapreads.tab", header = FALSE, stringsAsFactor = FALSE)

    read_counts$Run = gsub("data/00_preprocessing/kneaddata_bowtie/singlerun/","", read_counts$Run)
    read_counts = read_counts[!grepl("_2.fastq", read_counts$Run),]
    read_counts$Run = gsub("_1.fastq","", read_counts$Run)
    read_counts$Readcount = 2*read_counts$Readcount
    colnames(alned_ref) = c("Run", "Assembly")
    prop_alned_ref = merge(alned_ref, read_counts, by = "Run")

    colnames(alned_mag_sr)=c("Run", "MAGs")
    alned_mag_sr$`Assembly Approach`="Single-sample"

    # Combine alignment proportion tables from reference genomes and binned MAGs. 
    prop_alned_all = merge(prop_alned_ref, alned_mag_sr, by = "Run")
    prop_alned_all$percmags = prop_alned_all$MAGs/prop_alned_all$Readcount*100
    prop_alned_all$percassemb = prop_alned_all$Assembly/prop_alned_all$Readcount*100
    write.csv(prop_alned_all, "prop_reads_aligned.csv", quote = FALSE, row.names = FALSE)

    # Make figure
    ggplot(prop_alned_all, aes(x = percmags, y = percassemb, color=`Assembly Approach`)) +
        geom_point() + theme_classic() +
        xlab("Reads mapping to MAGs (%)") + ylab("Reads mapping to assembly (%)") +
        scale_color_manual(breaks=c("Single-sample","Coassembly"), values=c("#BBBBBB","#4477AA")) +
        xlim(0,100) + ylim(0,100)
    ggsave("prop_reads_aligned.png", width = 10, height = 10)
}

# setwd("./")
# source("~/Dropbox/workspace/MAG_Snakemake_wf/scripts/plotting/plot_ss.R")

# checkm_df = graph_checkm()
# graph_cmseq()
# graph_alignment()
# graph_dnadiff(checkm_df)
# graph_gtdb()
