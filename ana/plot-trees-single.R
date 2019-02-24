#library(tidyverse)
library(ggtree)
library(ape)
library(tools)

args = commandArgs(trailingOnly=TRUE)
in.file <- args[1]
branch.num <- as.numeric(args[2])

make.tree <- function(d){
   nedge <- nrow(d)
   nleaf <- (nedge + 2)/2
   nnode <- nleaf - 1

   mytree <- list()
   mytree$edge <- as.matrix(d[,c(1,2)])
   mytree$Nnode <- as.integer(nnode)
   mytree$tip.label <- paste(1:nleaf)
   mytree$edge.length <- d[,3]
   class(mytree) <- "phylo"
   checkValidPhylo(mytree)
   return(mytree)
}

plot.tree <- function(tree){
   p <- ggtree(tree) #+ geom_rootedge() 
   p <- p + geom_tiplab()
   p <- p + geom_text2(aes(subset=!isTip,label=node), hjust=-.3)
   edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge), edge_len=tree$edge.length)
   colnames(edge)=c("parent", "node", "edge_num", "edge_len")
   p <- p %<+% edge + geom_text(aes(x=branch, label=edge_len), nudge_y=0.1)
   print(p)
}

#dir <- "../sim-data/"
#dir <- data.dir
#files <- list.files(dir,"^sim\\-data\\-\\d+\\-tree")
files <- c(in.file)
for(f in files){
    cat("running on:", f, "\n")
    dd <- read.table(f,header=T)
    stub <- file_path_sans_ext(f)
    fout <- paste("plot-",stub,".pdf",sep="")
    cat("\n\nrunning over", stub, fout, "\n", sep="\t")
    
    #dd$start <- dd$start + 1
    #dd$end <- dd$end + 1

    if(branch.num==0){
        dd$length <- as.numeric(dd$length)
	dd <- dd[,c(1,2,3)]
    }
    if(branch.num==1){
	dd$nmut <- as.numeric(dd$nmut)
	dd <- dd[,c(1,2,5)]
    }

    mytree <- make.tree(dd)
    pdf(fout)
    plot.tree(mytree)
    dev.off()
}
