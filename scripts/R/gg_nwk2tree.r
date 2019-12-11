library("ape")
library("ggtree")
library("tidyverse")



# Takes one argument: filename: the filename of the txt file that contains a newick tree that we want a pdf made from.

args = commandArgs(TRUE)
###
fileName = args[length(args)]
fileName = "tree.nwk"
fileNameTrunc = tools::file_path_sans_ext(fileName)


file = readChar(fileName, file.info(fileName)$size)


tree <- read.tree(text = file)



length = max(tree$Nnode /18 * 210, 210)
p = ggtree(tree) + 
    geom_tiplab() + 
    geom_treescale(x = 0, width = signif(sum(tree$edge.length)*1, 1)) +
    labs(title = tools::file_path_sans_ext(fileName),
         subtitle = 'jModelTest2 ML')
p
ggsave(paste0(fileNameTrunc, ".pdf"), p, device = "pdf", width = length, height = length, units = "mm")
    

    
    


