library("ape")




# Takes 2 arguments
# 1: the input file, must be a newick formatted file
# 2: the output file, relative to the input file, must be <something>.pdf

args = commandArgs(TRUE)
fileName = args[length(args)-1]
fileout = args[length(args)]

print('these are the arguments:')
print(fileName)
print(fileout)



file = readChar(fileName, file.info(fileName)$size)


tree <- read.tree(text = file)

pdf(paste0("tree.pdf"), height = 10)#,width=6,height=4,paper='special')
plot(tree, cex = 0.5, label.offset = 0) #0.002
axisPhylo()
title(paste(tools::file_path_sans_ext(fileout)))


dev.off()

