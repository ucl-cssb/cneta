#library(tidyverse)
library(ggtree)
library(ape)
library(tools)

args = commandArgs(trailingOnly=TRUE)
data.file <- args[1]
annot.file <- args[2]
out.file <- args[3]
branch.num <- as.numeric(args[4])


make.tree <- function(d, labels){
   nedge <- nrow(d)
   nleaf <- (nedge + 2)/2
   nnode <- nleaf - 1

   mytree <- list()
   mytree$edge <- as.matrix(d[,c(1,2)])
   mytree$Nnode <- as.integer(nnode)
   #mytree$tip.label <- paste(1:nleaf)
   mytree$tip.label <- labels
   mytree$edge.length <- d[,3]
   class(mytree) <- "phylo"
   checkValidPhylo(mytree)
   return(mytree)
}

plot.tree <- function(tree){
   p <- ggtree(tree, size=0.5, linetype=1) #+ geom_rootedge() 
   p <- p + geom_tiplab(align=TRUE) + theme_tree2() + xlim(NA, 130)
   p <- p + geom_text2(aes(subset=!isTip,label=node), hjust=-.3)
   edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge), edge_len=tree$edge.length)
   colnames(edge)=c("parent", "node", "edge_num", "edge_len")
   p <- p %<+% edge + geom_text(aes(x=branch, label=edge_len), nudge_y=0.1)
   print(p)
}


dd <- read.table(data.file,header=T)
da <- read.table(annot.file,header=T)
da$sample <- as.character(da$sample)


if(branch.num==0){
   dd$length <- as.numeric(dd$length)
   dd[ which(dd$length < 1e-3), ]$length <- 0

   dd <- dd[,c(1,2,3)]
}
if(branch.num==1){
   dd$nmut <- as.numeric(dd$nmut)
   dd <- dd[,c(1,2,5)]
}

mytree <- make.tree(dd, as.character(da$sample))
pdf(out.file)
plot.tree(mytree)
dev.off()

