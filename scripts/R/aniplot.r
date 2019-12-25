library(tidyverse)
library(reshape2)

# Usage: Rscript aniplot.r <path_to_tsv> <path_to_newick_tree>


# Import arguments
args = commandArgs(TRUE)
file_name = args[length(args)-1] # second last argument
tree_file = args[length(args)] # last argument


# extract order from tree
tree = readChar(tree_file, file.info(tree_file)$size)

newick2names = function(newick) {
  splitted = newick %>% str_split(",")
  
  
  rv = c()
  for (i in splitted) {
    
    new_i = str_replace(i, ":.+$", "") %>% str_replace("\n", "") %>% str_replace("\\(+", "") %>% str_replace("\\)+", "")
    print(new_i)
    rv = c(rv, new_i)
  }
  rv
}

tree_as_vector = newick2names(tree)

# get truncated file name
file_name_trunc = tools::file_path_sans_ext(file_name)

# read tsv-file
ani <- read_delim(file_name, "\t", escape_double = FALSE, trim_ws = TRUE, na = "-")

# melt for plotting
melted_ani <- melt(ani, na.rm = T)
melted_ani_2 = tibble(X1 = melted_ani$variable,
                      variable = melted_ani$X1,
                      value = melted_ani$value)


font_size = 5

p = ggplot(data = melted_ani, aes(x=X1, y=variable, fill=value)) + 
    geom_tile() +
    #geom_tile(data = rbind(melted_ani_2), aes(x=X1, y = variable, fill = value))+
    geom_tile(aes(x = variable, y = X1, fill = value)) +
    #xlim(rev(levels(melted_ani$variable))) + 
    #ylim(rev(levels(melted_ani$variable))) +
    xlim(tree_as_vector) + 
    ylim(tree_as_vector) +
    labs(title = file_name_trunc,
         subtitle = "Average Nucleotide Identity of core gene alignment (sanger-panito)",
         x = "samples",
         y = "samples") +
    scale_fill_gradient(name = "ANI [%]") + 
    #geom_text(aes(label = round(value, 1)), size = 2, alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 90, size = font_size, hjust = TRUE),
          axis.text.y = element_text(size = font_size))

#p

# Save to disk next to tsv-file.
ggsave(paste0(file_name_trunc, ".pdf"), p, device = "pdf", width = 210, height = 210, units = "mm")





print(tree_as_vector)

