library(plyr) # Must be loaded before to ensure correct overloadings
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(tidyr)
library(tools)


# ID_TOOL_matches.csv: Table [binner,first,second,dist] of best matches between bins made using different workflows.
# ID_TOOL_distances.csv: CSV of reciprocal distances between bins made using different workflows.
# ID_TOOL_sumstats.csv: CSV of the number of reciprocal, non-reciprocal, and no-match bins compared to the entire bin infere

#################
# Single metric #
#################

graph_distances <- function(pt_id, tool) {
    match_df = read.csv(paste(pt_id, tool, 'matches.csv', sep = "_"), header = TRUE, stringsAsFactor = FALSE)
    dist_df = transpose(read.csv(paste(pt_id, '_', tool, '_distances.csv', sep = ""), header = FALSE, stringsAsFactor = FALSE))
    dist_df$reciprocal = "All"
    colnames(dist_df)=c("dist", "reciprocal")
    ggplot(dist_df, aes(x = dist, color = reciprocal, fill = reciprocal)) +
        geom_histogram(alpha = 0.3) +
        geom_histogram(data = match_df, alpha = 0.3) +
        xlab(toTitleCase(tool)) + ylab("Count") +
        guides(color = FALSE) + labs(fill="Reciprocal\nbest-match?")
    ggsave(paste(pt_id, tool, "distances.png", sep = "_"), width = 15, height = 10)
}


process_df <- function(id, tool) {
    df = read.csv(paste(id, tool, 'matches.csv', sep = "_"), header = TRUE, stringsAsFactor = FALSE)
    df[ , tool] <- NA
    for(i in 1:nrow(df)) {
        binner = df[i,1]
        if (binner == 'B') {
            other = 'H'
        } else { other = 'B' }
        df[i,2] <- paste(binner, df[i,2], sep = '_')
        df[i,3] <- paste(other, df[i,3], sep = '_')
        df[i,5] <- as.logical(df[i,5])
    }
    names(df) = c('binner', 'bin', tool, paste(tool, 'dist', sep = "_"), paste(tool, 'recip', sep = "_"))
    return(df[,2:5])
}


###################
### Two metrics ###
###################


graph_agreement <- function(pt_id, tool_one, tool_two) {
    df_one = process_df(pt_id, tool_one)
    df_two = process_df(pt_id, tool_two)
    df_compare = merge(df_one, df_two, by = 'bin') # Intersect df on the bin

    recip_one = paste(tool_one, 'recip', sep = "_")
    recip_two = paste(tool_two, 'recip', sep = "_")
    df_compare$r_best_match = ifelse(df_compare[recip_one] == TRUE & df_compare[recip_two] == TRUE,TRUE, FALSE)
    df_compare$metric_agree = ifelse(df_compare[tool_one]==df_compare[tool_two],TRUE,FALSE)
    write.csv(df_compare, paste(pt_id, "agreements.csv", sep = '_'), quote = FALSE, row.names = FALSE)

    ct_compare = as.data.frame(table(df_compare[,c("r_best_match", "metric_agree")]))
    names(ct_compare) = c("Reciprocal Best Matches", "Metric Agreement","Count")
    ggplot(ct_compare, aes(x=`Reciprocal Best Matches`, y=`Metric Agreement`, fill=Count)) + 
        geom_tile() + coord_equal() +
        geom_text(aes(label=Count), color="black") + # Print count in square
        theme(legend.position="none")
    ggsave(paste(pt_id, "_agreements.png", sep = ''),width=15)
}


# setwd("./")
# source("~/Dropbox/workspace/MAG_Snakemake_wf/scripts/plotting/plot_m2m.R")

# ids = c('98', '99', '101', '102')
# metrics = c('fastani', 'mash')
# for(i in 1:length(ids)) {
#     graph_distances(ids[i], metrics[1])
#     graph_distances(ids[i], metrics[2])
#     graph_agreement(ids[i], metrics[1], metrics[2])
# }
