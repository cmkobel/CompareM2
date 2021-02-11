
library(tidyverse)
library(ggtree)
#library(tidytree)

args = commandArgs(trailingOnly=TRUE)
paste("arguments given:", args)


newick_file = args[1]
mldist_file = args[2]
title = args[3]
#setwd("C:/Users/carkob/Desktop/bootstraptree/iq3/"); newick_file = "simple.fasta.treefile"; title = "Some title"; mldist_file = "simple.fasta.mldist"

cat("newick_file:", newick_file, "\n")
cat("title:", title)
cat("mldist_file:", mldist_file, "\n")

tree <- read.tree(newick_file)


height = max(7, length(tree$tip.label)/1.2) #minimum 7
height = min(height, 50) # maximum 50

width = max(6, height/1.5)


p1 = ggtree(tree) +
    geom_tree() +
    theme_tree2() +
    geom_tiplab(align=FALSE, color = "grey30") +
    
    geom_nodelab(color = "grey40", hjust = 0) +
    labs(title = title,
         subtitle = "concatenated core genome alignment maximum likelihood phylogeny (unrooted)",
         caption = "internal nodes: percentage bootstrap support of 1000 replicates.\nIQ-TREE2 (B. Q. Minh et al. 2020)")
new_max = ggplot_build(p1)$layout$panel_scales_x[[1]]$range$range[2] 

p2 = p1 + xlim(0, new_max*1.12)
p2 +



ggsave(paste0("iqtree_ml_bootstrap.pdf"), plot = p2, width = width, height = height, limitsize = F)



# Histogram of distances
mldist = simple_fasta <- read_table2(mldist_file, 
    col_names = FALSE, skip = 1)


#trim 
mldist = mldist[1:(dim(mldist)[1]+1)] %>% column_to_rownames("X1")
mldist[lower.tri(mldist)] <- NA
mldist[!upper.tri(mldist)] <- NA
#mldist = mldist[-1]

mldist = rownames_to_column(mldist, "from")


names(mldist) = c("from", mldist$from)


tmp = mldist %>% pivot_longer(-from, names_to = "to", values_to = "distance") %>% 
    filter(from != to) %>%
    drop_na() %>% mutate(pairs = paste0(from, '\n', to)) %>% arrange(distance)
    
# auto hist
tmp %>% ggplot(aes(distance, fill = pairs)) + 
    geom_histogram() +
    # stat_bin(bins=30, geom="text", colour="black",
    #        aes(label=pairs, group = pairs),
    #        position=position_stack(vjust=0.5)
    #        ) +
    #stat_bin(bins = 30, geom='text', color='white',
    #       position=position_stack(vjust = 0.5)) +
    #geom_text(stat = "count", aes(label = pairs, y = stat(count)), position=position_stack(0.5), angle = 90)+
    #geom_text(aes(x = distance, y = stat(count), label = pairs))
    labs(title = "Histogram of distances", 
         x = "maximum likelihood distance") +
ggsave(paste0(title, "_histogram.pdf"), limitsize = F)


# manual hist
#space = (max(tmp$distance) - min(tmp$distance))/30
a = tmp %>% mutate(group = cut_width(distance, (max(distance)-min(distance))/30)) %>% 
    group_by(group) %>% 
    mutate(distance = mean(distance), count = row_number())

a %>% summarize(dist = mean(distance), n = length(distance)) %>% 
    ggplot(aes(dist, n)) +
    geom_col() +
    geom_text(data = a,
                 mapping = aes(label = pairs,
                               x = distance,
                               y = count-0.95),
              angle = 90, hjust = 0, vjust = 0.5, size = 1.2, color = "white") +
    labs(title = "Histogram of distances", x = "maximum likelihood distance", y = "count") +


    # geom_segment(data = a,
    #              mapping = aes(x = distance, xend = distance,
    #                            y = count-1, yend = count),
    #              color = "red", width = space)

    
ggsave(paste0(title, "_histogram2.pdf"), limitsize = F)
    


tmp %>% ggplot(aes(x = pairs, y = distance)) + 
    geom_col() +
    scale_x_discrete(limits = tmp$pairs) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = "Sorted pairwise distances",
         y = "maximum likelihood distance") +
ggsave(paste0(title, "_distances.pdf"), limitsize = F)

