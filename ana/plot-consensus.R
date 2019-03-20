#!/usr/bin/env Rscript

library(ggtree)
library(ape)


args = commandArgs(trailingOnly=TRUE)
dir.bstrap = args[1]
file.tree = args[2]
file.out = args[3]
file.sample = NA
if(length(args)>3) file.sample = args[4]

cat(dir.bstrap,file.tree,file.out,"\n")

# dir.bstrap = "D:/Gdrive/git/sveta/test/bootstrap"
# file.tree="D:/Gdrive/git/sveta/test/results-maxL-tree-sim1.txt"
# file.out="D:/Gdrive/git/sveta/test/results-maxL-tree-sim1-bootstap.pdf"

make.tree <- function(d, labels){
   nedge <- nrow(d)
   nleaf <- (nedge + 2)/2
   nnode <- nleaf - 1

   mytree <- list()
   mytree$edge <- as.matrix(d[,c(1,2)])
   mytree$Nnode <- as.integer(nnode)
   if(is.na(labels)){
     mytree$tip.label <- paste(1:nleaf)
   }
   else{
     mytree$tip.label <- labels
   }
   mytree$edge.length <- d[,3]
   class(mytree) <- "phylo"
   checkValidPhylo(mytree)
   return(mytree)
}


plot.tree <- function(tree){
  p <- ggtree(tree) #+ geom_rootedge()
  # support <- character(length(tree$node.label))
  # #The following three lines define your labeling scheme.
  # support[tree$node.label >= 95] <- "red"
  # support[tree$node.label < 95 & tree$node.label >= 70] <- "pink"
  # support[tree$node.label < 70] <- "blue"
  p <- p + geom_tiplab() #+ xlim(c(0,60))
  p <- p + geom_text2(aes(subset=!isTip, label=label, hjust=-.3, color="red"))
  edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge), edge_len=tree$edge.length)
  colnames(edge)=c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x=branch, label=edge_len), nudge_y=0.1)
  print(p)

}


dt = read.table(file.tree,header=T)
labels = NA
if(!is.na(file.sample)){
  da = read.table(file.sample,header=T,stringsAsFactors = F)
  labels = da$sample
}
tree = make.tree(dt, labels)

btrees=list()
files = list.files(path = dir.bstrap, pattern = "*.txt", recursive=F)
for (i in 1:length(files)){
  fname = file.path(dir.bstrap, files[i])
  dt = read.table(fname,header=T)
  btree = make.tree(dt, labels)
  btrees[[i]] = btree
}
# mtree = consensus(btrees, p=0.5)

# get proportions of each clade:
clad <- prop.clades(tree, btrees, rooted = TRUE)
clad[is.na(clad)] = 0
tree$node.label = clad * 100 / length(files)

# pID="1"
# title = paste0("Patient ", pID, " - ML tree with bootstrap support")
# plot(tree, main = title)
#plot.tree(tree)
plot.tree(tree)
ggsave(file.out, width=11.69, height=8.27, units="in", limitsize = FALSE)
# plot(mtree)
