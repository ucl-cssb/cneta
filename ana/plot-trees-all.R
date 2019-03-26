# library(tidyverse)
library(ggtree)
library(ape)
library(tools)

if(!require(optparse)){
    install.packages("optparse")
    library(optparse)
}

# This script is used to plot phylogenetic trees:
# 1) plotting a single tree file with tips as 1 to n or with provided labels
# 2) plotting a single tree file with bootstrapping support
# 3) plotting all tree files in a directory

make.tree <- function(d, labels = NA) {
  nedge <- nrow(d)
  nleaf <- (nedge + 2)/2
  nnode <- nleaf - 1

  mytree <- list()
  mytree$edge <- as.matrix(d[, c(1, 2)])
  mytree$Nnode <- as.integer(nnode)

  if(is.na(labels)){
     mytree$tip.label <- paste(1:nleaf)
  }
  else{
    mytree$tip.label <- labels
  }

  mytree$edge.length <- d[, 3]
  class(mytree) <- "phylo"
  checkValidPhylo(mytree)
  return(mytree)
}


plot.tree <- function(tree) {
  p <- ggtree(tree)  #+ geom_rootedge()
  p <- p + geom_tiplab()
  p <- p + geom_text2(aes(subset = !isTip, label = node), hjust = -0.3)
  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge) = c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)
  print(p)
}

# Plot tree with xlim specified to show full tip labels
plot.tree.xlim <- function(tree) {
  p <- ggtree(tree, size = 0.5, linetype = 1)  #+ geom_rootedge()
  # Add margin to show full name of labels  if (is.na(tree.max))
  tree.max = max(node.depth.edgelength(tree)) + 20
  p <- p + geom_tiplab(align = TRUE) + theme_tree2() + xlim(NA, tree.max)
  # p <- p + geom_text2(aes(subset=!isTip,label = node), hjust=-.3)
  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge) = c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)
  print(p)
}


plot.tree.bootstrap <- function(tree){
  p <- ggtree(tree) #+ geom_rootedge()
  # support <- character(length(tree$node.label))
  # #The following three lines define your labeling scheme.
  # support[tree$node.label >= 95] <- "red"
  # support[tree$node.label < 95 & tree$node.label >= 70] <- "pink"
  # support[tree$node.label < 70] <- "blue"
  tree.max= max(node.depth.edgelength(tree)) + 20
  p <- p + geom_tiplab(align = TRUE) + theme_tree2() + xlim(NA, tree.max)
  p <- p + geom_text2(aes(subset=!isTip, label=label, hjust=-.3, color="red"))
  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge)=c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)
  print(p)
}

print.tree <- function(tree_file, tree_style, out_file="", branch_num = 0, bstrap_dir="", pattern="", labels = NA) {
  dd <- read.table(tree_file, header = T)
  dir <- dirname(tree_file)
  stub <- file_path_sans_ext(basename(tree_file))
  if(out_file!=""){
    fout <- out_file
  }
  else{
    fout <- file.path(dir, paste("plot-",stub,".pdf",sep=""))
  }
  cat("\n\nrunning over", stub, fout, "\n", sep = "\t")

  # dd$start <- dd$start + 1 dd$end <- dd$end + 1
  if (branch_num == 0) {
    dd$length <- as.numeric(dd$length)
    dd[which(dd$length < 1e-3), ]$length <- 0
    dd <- dd[, c(1, 2, 3)]
  }
  if (branch_num == 1) {
    dd$nmut <- as.numeric(dd$nmut)
    dd <- dd[, c(1, 2, 5)]
  }

  mytree <- make.tree(dd, labels)

  pdf(fout)
  if(tree_style=="simple") plot.tree(mytree)
  if(tree_style=="xlim") plot.tree.xlim(mytree)
  if(tree_style=="bootstrap"){  # bootstrap
    btrees = list()
    cat("Patterns to match bootstrapping trees: ", pattern, "\n")
    files = list.files(path = bstrap_dir, pattern = glob2rx(pattern), recursive = F)
    for (i in 1:length(files)){
      fname = file.path(bstrap_dir, files[i])
      dt = read.table(fname,header = T)
      btree = make.tree(dt, labels)
      btrees[[i]] = btree
    }
    clad <- prop.clades(mytree, btrees, rooted = TRUE)
    clad[is.na(clad)] = 0
    mytree$node.label = clad * 100 / length(files)
    plot.tree.bootstrap(mytree)
  }
  dev.off()
  # ggsave(file.out, width = 11.69, height = 8.27, units="in", limitsize = FALSE)
}

get.labels <- function(annot_file){
  labels = NA
  if(annot_file!=""){
    da <- read.table(annot_file,header = T,stringsAsFactors = F)
    labels = as.character(da$sample)
  }
  return(labels)
}


option_list = list(
  make_option(c("-f", "--tree_file"), type="character", default="",
              help="dataset file name [default=%default]", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default="",
              help="The name of output file [default=%default]", metavar="character"),
  make_option(c("-a", "--annot_file"), type="character", default="",
              help="The file containing the labels of tip nodes [default=%default]", metavar="character"),
  make_option(c("-d", "--tree_dir"), type="character", default="",
              help="The directory containing all the tree files to plot [default=%default]", metavar="character"),
  make_option(c("-s", "--bstrap_dir"), type="character", default="",
              help="The directory containing all the bootstrapping tree files [default=%default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="",
              help="The naming pattern of tree files [default=%default]", metavar="character"),
	make_option(c("-b", "--branch_num"), type="integer", default = 0,
              help="The type of values on branches (0: branch length, 1: mutation number) [default=%default]", metavar="integer"),
  make_option(c("-t", "--plot_type"), type="character", default="single",
              help="The type of plot, including: all (plotting all tree files in a directory), single (plotting a single tree file), and bootstrap (plotting a single tree file with bootstrapping support) [default=%default]", metavar="character"),
  make_option(c("-l", "--tree_style"), type="character", default="simple",
              help="The style of tree plot, including: simple (a simple tree with tip labels and branch lengths), xlim (adding xlim to the tree), and bootstrap (adding bootstrap support value to the tree) [default=%default]", metavar="character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

tree_file = opt$tree_file
out_file = opt$out_file
plot_type = opt$plot_type
tree_dir = opt$tree_dir
pattern = opt$pattern
bstrap_dir = opt$bstrap_dir
tree_style = opt$tree_style
annot_file = opt$annot_file
branch_num = opt$branch_num

# cat("Parameters used here:\n")
# cat("tree_file:", tree_file, "\n")
# cat("plot_type:", plot_type, "\n")
# cat("tree_dir:", tree_dir, "\n")
# cat("branch_num:", branch_num, "\n")
# cat("pattern:", pattern, "\n")

if (plot_type == "all"){
  # dir <- '../sim-data/'
  dir <- tree_dir
  cat(paste0("Plotting all trees in directory ", dir, "\n"))
  if(pattern == ""){
    files <- list.files(dir, "^sim\\-data\\-\\d+\\-tree")
  }
  else{
    files <- list.files(dir, pattern = glob2rx(pattern))
  }
  #print(files)

  for (f in files) {
    cat("running on:", f, "\n")
    fname = file.path(dir, f)
    print.tree(fname, tree_style, branch_num = branch_num)
  }
}

if (plot_type == "single"){
  cat(paste0("Plotting the tree in ", tree_file))
  labels = get.labels(annot_file)
  print.tree(tree_file, tree_style, out_file = out_file, branch_num = branch_num, labels = labels)
}

if (plot_type == "bootstrap"){
  cat(paste0("Plotting bootstrap values for the tree in ", tree_file))
  tree_style = "bootstrap"
  labels = get.labels(annot_file)
  print.tree(tree_file, tree_style, out_file = out_file, branch_num = branch_num, bstrap_dir = bstrap_dir, pattern = pattern, labels = labels)
}
