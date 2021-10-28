library(plyr) # Must be loaded before to ensure correct overloadings
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(tidyr)

# metaQUAST_report.tsv: Report of reference-free statistics from metaQUAST.
# binning_stats.csv: Basic binning statistics of single-sample and coassembly refined bins. Comes with binning_stats.lst, which indicates which bins came from which runs. 
# Figure 3a,b
# sr_summ_cmseq_all.csv: CMSeq (strain heterogeneity) summary statistics of single-sample MAGs. 
# coas_summ_cmseq_all.csv: CMSeq (strain heterogeneity) summary statistics of coassembled MAGs. 
# # Figure 3c
# sr_checkm_metrics.csv: Completeness, contamination, and strain heterogeneity of each single sample-assembled refined MAG bin.
# coas_checkm_metrics.csv: Completeness, contamination, and strain heterogeneity of each coassembled refinedMAG bin.
# Figure 3d,e
# dnadiff_summary.tsv: Comparison between each refined MAG bin (single-sample and coassembly) and the most similar RefSeq reference genome. Includes total bases in reference genome, the number of aligned bases, and the average identity. 
# Figure 3f
# gtdbtk.bac120.summary.tsv: GTDB classifications of each classifiable refined MAG bin, with the ANI to the most similar RefSeq reference genome. 
# Figure 4a
# readcounts.tsv: Number of reads in each sample's FastQ after Kneaddata cleaning. 
# flagstat_sum.txt: Proportion of reads in each sample's FastQ that are either uniquely mapped or chimeric over all mapped reads to the assembled scaffolds (the other category is multi-mapped).
# sr_catalogue_mapreads.tab: Proportion of reads in each sample's FastQ that are either uniquely mapped or chimeric over all mapped reads to the refined MAG bins (the other category is multi-mapped).
# coas_catalogue_mapreads.tab: sr_catalogue_mapreads.tab: Proportion of reads in the FastQs in each coassembly set that are either uniquely mapped or chimeric over all mapped reads to the refined MAG bins (the other category is multi-mapped).

# Figure 3a,b: Box-and-whisker plot comparing completeness and contamination.
graph_checkm <- function(){
    # Load dataframes and format column names/variable types
    checkm_sr = read.csv("sr_checkm_metrics.csv", header = TRUE, stringsAsFactor = FALSE)
    checkm_co = read.csv("coas_checkm_metrics.csv", header = TRUE, stringsAsFactor = FALSE)
    checkm_sr$status="Single-sample"
    checkm_co$status="Coassembly"

    # Combine single-sample and coassembly tables
    checkm = rbind.data.frame(checkm_sr, checkm_co)
    checkm$status = factor(checkm$status, levels = c("Single-sample", "Coassembly"))

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
    cmseq_sr = read.csv("sr_summ_cmseq_all.csv", header = TRUE, stringsAsFactor = FALSE)
    cmseq_co = read.csv("coas_summ_cmseq_all.csv", header = TRUE, stringsAsFactor = FALSE)
    cmseq_sr$status="Single-sample"
    cmseq_co$status="Coassembly"
    colnames(cmseq_sr)=c("path", "strain_het")
    colnames(cmseq_co)=c("path", "strain_het")
    cmseq_sr$strain_het=as.vector(cmseq_sr$strain_het))
    cmseq_co$strain_het=as.vector(cmseq_co$strain_het))
    cmseq_sr=cmseq_sr[complete.cases(cmseq_sr$strain_het),]
    cmseq_co=cmseq_co[complete.cases(cmseq_co$strain_het),]

    # Combine single-sample and coassembly tables
    cmseq = rbind.data.frame(cmseq_sr, cmseq_co)
    cmseq$status = factor(cmseq$status, levels=c("Single-sample", "Coassembly"))

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
          theme(axis.title.x = element_text(size=14) + axis.text.x = element_text(size=12))) + 
          theme(axis.title.y = element_blank() + axis.text.y = element_text(size=12))
    ggsave("gtdb_classification.png", tax_class, width = 15, height = 10)
}

# Figure 4a: Dot plot of reads aligned to the identified reference genomes vs. binned MAGs.
graph_alignment <- function() {
    # Load dataframes and format column names/variable types
    read_counts = read.table("readcounts.tsv", header = TRUE, stringsAsFactor = FALSE)
    alned_ref = read.table("flagstat_sum.txt", header = FALSE, stringsAsFactor = FALSE)
    alned_mag_sr = read.table("sr_catalogue_mapreads.tab", header = FALSE, stringsAsFactor = FALSE)
    alned_mag_co = read.table("coas_catalogue_mapreads.tab", header = FALSE, stringsAsFactor = FALSE)

    read_counts$Run = gsub("data/00_preprocessing/kneaddata_bowtie/singlerun/","", read_counts$Run)
    read_counts = read_counts[!grepl("_2.fastq", read_counts$Run),]
    read_counts$Run = gsub("_1.fastq","", read_counts$Run)
    read_counts$Readcount = 2*read_counts$Readcount
    colnames(alned_ref) = c("Run", "Assembly")
    prop_alned_ref = merge(alned_ref, read_counts, by = "Run")

    colnames(alned_mag_sr)=c("Run", "MAGs")
    alned_mag_sr$`Assembly Approach`="Single-sample"
    colnames(alned_mag_co)=c("Run", "MAGs")
    alned_mag_co$`Assembly Approach`="Coassembly"
    prop_alned_mag = rbind.data.frame(alned_mag_sr, alned_mag_co)

    # Combine alignment proportion tables from reference genomes and binned MAGs. 
    prop_alned_all = merge(prop_alned_ref, prop_alned_mag, by = "Run")
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

horizontal_dot <- function(df, cname) {
    p = ggplot(df, aes(x = value, y = variable, colour = eval(as.symbol(cname)))) + 
        geom_point(size = 5) + # theme_classic() + 
        theme(axis.title.x = element_blank(), axis.text.x = element_text(size = rel(1.2))) + 
        theme(axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.4))) + 
        guides(color=guide_legend(title="Approach")) + theme(legend.position="none")
    if (max(df$value) > 10000) {
        p = p + scale_x_continuous(n.breaks = 6, labels = scales::comma) # Increases the number of x-axis ticks
    } else {
        p = p + scale_x_continuous(n.breaks = 6)
    }
    return(p)
}

# Figures: De novo assembly summary statistics based on metaQUAST's reference-free report.
graph_quast <- function() {
    # Load dataframes and format column names/variable types    
    raw_df = read.csv("metaQUAST_report.tsv", header = FALSE, row.names = 1, stringsAsFactor = FALSE)
    parameter_names = rownames(raw_df)
    trn_df = transpose(raw_df)
    quast <- as.data.frame(lapply(trn_df, function(x) as.numeric(x)))
    colnames(trn_df) = parameter_names
    colnames(quast) = parameter_names
    quast$Assembly = trn_df$Assembly
    quast$Assembly = gsub("_scaffolds","", quast$Assembly)

    # Make average usage and average contig length columns
    quast$Avg_ctg_used = quast$Num_contigs / quast$Num_contigs_0_bp
    quast$Avg_len_used = quast$Total_length / quast$Total_length_0_bp
    quast$Avg_binned_contig = quast$Total_length / quast$Num_contigs
    quast$Avg_all_contig = quast$Total_length_0_bp / quast$Num_contigs_0_bp

    # Reshape dataframe such that the parameters are variables
    quast_mod = melt(quast, id.vars=c("Assembly"))
    quast_mod$value = as.numeric(quast_mod$value)

    # Plot 1a: Total (binnable) assembled length 
    plot_1a = horizontal_dot(quast_mod[quast_mod$variable == "Total_length",], "Assembly") + 
        scale_y_discrete(labels=c("Total_length" = "Total\nAssembly\nLen. (bp)"))

    # Plot 1b: Number of (binnable) contigs
    plot_1b = horizontal_dot(quast_mod[quast_mod$variable == "Num_contigs",], "Assembly") +
        scale_y_discrete(labels=c("Num_contigs" = "Number of\nContigs"))

    # Plot 1c,d: Proportion of binnable vs. unused contigs and lengths
    plot_1c = horizontal_dot(quast_mod[quast_mod$variable == "Avg_len_used",], "Assembly") + 
        scale_y_discrete(labels=c("Avg_len_used" = "Proportion of\nBinnable\nAssembly"))
    plot_1d = horizontal_dot(quast_mod[quast_mod$variable == "Avg_ctg_used",], "Assembly") + 
        scale_y_discrete(labels=c("Avg_ctg_used" = "Proportion of\nBinnable\nContigs")) + 
        theme(legend.position="bottom", legend.title=element_text(size=rel(1.2)), legend.text=element_text(size=rel(1.0)))

    # Plot 2: Largest binnable, average binnable, and average all contig length
    plot_2a = horizontal_dot(quast_mod[quast_mod$variable == "Largest_contig",], "Assembly") +
        scale_y_discrete(labels=c("Largest_contig" = "Longest\nBinnable\nContig (bp)"))
    plot_2b = horizontal_dot(quast_mod[quast_mod$variable == "Avg_binned_contig",], "Assembly") +
        scale_y_discrete(labels=c("Avg_binned_contig" = "Average\nBinnable\nContig (bp)"))
    plot_2c = horizontal_dot(quast_mod[quast_mod$variable == "Avg_all_contig",], "Assembly") +
        scale_y_discrete(labels=c("Avg_all_contig" = "Average of All\nContigs (bp)")) +
        theme(legend.position="bottom", legend.title=element_text(size=rel(1.2)), legend.text=element_text(size=rel(1.0)))

    # Plot 3: NX and LX
    plot_3a = horizontal_dot(quast_mod[quast_mod$variable == "N50" | quast_mod$variable == "N75",], "Assembly")
    plot_3b = horizontal_dot(quast_mod[quast_mod$variable == "L50" | quast_mod$variable == "L75",], "Assembly") +
        theme(legend.position="bottom", legend.title=element_text(size=rel(1.2)), legend.text=element_text(size=rel(1.0)))

    # Make multi-panel figures
    png("assembly_sumstats.png", width = 900, height = 400)
    grid.draw(rbind(ggplotGrob(plot_1a), ggplotGrob(plot_1c), ggplotGrob(plot_1b), ggplotGrob(plot_1d)))
    dev.off()
    png("contig_sumstats.png", width = 900, height = 300)
    grid.draw(rbind(ggplotGrob(plot_2a), ggplotGrob(plot_2b), ggplotGrob(plot_2c)))
    dev.off()
    png("nx_lx_sumstats", width = 900, height = 200)
    grid.draw(rbind(ggplotGrob(plot_3a), ggplotGrob(plot_3b)))
    dev.off()
}

# Figures: MAG bin N50s and sizes and the binning algorithms they were generated from.
graph_binning <- function() {

    # Load bin count dataframe and make a vector of sample names
    count_df = read.table("binning_stats.lst")
    count_df$V2 = gsub("data/02_binning/", "", count_df$V2)
    count_df$V2 = gsub("/metawrap_bin_refinement/metawrap_50_10_bins.stats", "", count_df$V2)
    count_df$V2 = gsub("coassembly/", "", count_df$V2)
    count_df$V2 = gsub("singlerun/", "", count_df$V2)
    count_clean_df = head(count_df, -1)
    samples = c()
    for (i in (1:nrow(count_clean_df))) {
        samples = append(samples, replicate(count_clean_df$V1[i] - 1, count_clean_df$V2[i]))
    }

    # Load statistics dataframe and add the sample name vector. Reshape and set variable types
    stats_df = read.csv("binning_stats.csv", header = TRUE)
    stats_df$sample = samples
    stats_mod_df = melt(stats_df, id.vars=c("sample"))
    stats_mod_df$sample = as.factor(stats_mod_df$sample)
    stats_mod_df$value = as.numeric(stats_mod_df$value)

    # Make N50 and bin size multi-panel figure
    plot_a = horizontal_dot(stats_mod_df[stats_mod_df$variable == "N50",], "sample")
    plot_b = horizontal_dot(stats_mod_df[stats_mod_df$variable == "size",], "sample") + 
        scale_y_discrete(labels=c("size" = "Bin Size\n(bp)"))

    png("bin_n50_size.png", width = 900, height = 200)
    grid.draw(rbind(ggplotGrob(plot_a), ggplotGrob(plot_b))
    dev.off()

    # Make standalone vs. refiner table from the bin column in the statistics dataframe, then figure 
    binners_ref = stats_df$binner
    binners_ref = gsub("bins", "", binners_ref)
    bin_step_df = as.data.frame(nchar(binners_ref) > 1) # Bins that came from single binner or refinement
    colnames(bin_step_df) = c("step")
    bin_step_df$step = revalue(as.character(bin_step$step), c("FALSE"="Standalone Binner", "TRUE"="Bin Refiner"))

    ggplot(bin_step_df, aes(x = step, fill = step)) + geom_bar() +
        ylab("Number of Bins") + xlab("Binning Step") +
        guides(color=guide_legend(title="Approach")) + theme(legend.position="none")
    ggsave("standalone_refinement_bins.png", width = 10, height = 10)

    # Make table of standalone binners, then figure 
    binners_std = strsplit(paste(binners_ref, collapse = ""), split = "")[[1]]  # Paste puts elements in a vector together into a string, and strsplit pulls them apart into chars
    bin_std_df = as.data.frame(binners_std)
    bin_std_df$binners_std = revalue(as.character(bin_std_df$binners_std), c("A"="MetaBAT(2)", "B"="MaxBin2", "C"="CONCOCT"))

    ggplot(bin_std_df, aes(x = binners_std, fill = binners_std)) + geom_bar() +
        ylab("Number of Bins") + xlab("Binning Algorithm") +
        guides(color=guide_legend(title="Approach")) + theme(legend.position="none")
    ggsave("standalone_binner.png", width = 10, height = 10)
}

# setwd("./")
# source("~/Dropbox/workspace/MAG_Snakemake_wf/scripts/plotting/plot_all.R")

# checkm_df = graph_checkm()
# graph_cmseq()
# graph_alignment()
# graph_dnadiff(checkm_df)
# graph_gtdb()
# graph_quast()
# graph_binning()
