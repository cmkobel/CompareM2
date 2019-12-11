library("ape")




# Takes one argument: filename: the filename of the txt file that contains a newick tree that we want a pdf made from.

args = commandArgs(TRUE)
fileName = args[length(args)]



file = readChar(fileName, file.info(fileName)$size)


tree <- read.tree(text = file)

pdf(paste0(fileName, ".pdf"), height = 10)#,width=6,height=4,paper='special')
plot(tree, cex = 0.5, label.offset = 0) #0.002
axisPhylo()
title(paste(tools::file_path_sans_ext(fileName)))


dev.off()

