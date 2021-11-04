#!/usr/bin/env Rscript

suppressMessages(library(ggtree))
suppressMessages(library(treeio))
suppressMessages(library(ape))
suppressMessages(library(tidyverse))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(optparse))


# This script is used to plot phylogenetic trees:
# 1) plotting a single tree file with tips as 1 to n or with provided labels
# 2) plotting a single tree file with bootstrapping support
# 3) plotting all tree files in a directory

MIN_BLEN = 1e-3

# d is a data frame with 3 columns: start, end, branch length
# Leaves must be encoded from 1 to nleaf
make.tree <- function(d, labels = c(), digit = 2) {
  nedge <- nrow(d)
  nleaf <- (nedge + 2) / 2
  nnode <- nleaf - 1

  mytree <- list()
  mytree$edge <- as.matrix(d[, c("start", "end")])
  mytree$Nnode <- as.integer(nnode)

  if(length(labels) < nleaf){
    mytree$tip.label <- paste(1:nleaf)
  }else{
    mytree$tip.label <- labels
  }

  # control precision of branch length
  mytree$edge.length <- round(d[, 3], digit)
  class(mytree) <- "phylo"
  # checkValidPhylo(mytree)

  return(mytree)
}


get.outfile.name <- function(tree_file){
  dir <- dirname(tree_file)
  stub <- file_path_sans_ext(basename(tree_file))
  
  if(branch_num == 1){
    mfix = paste0(stub, "-nmut")
  }else{
    mfix = stub
  }
  
  fout <- file.path(dir, paste("plot-", mfix, ".pdf", sep=""))

  cat("\n\nrunning over", stub, fout, "\n", sep = "\t")
  
  return(fout)
}

# Prepare the tree for plotting
get.tree <- function(tree_file, branch_num = 0, labels = c(), scale_factor = 1){
  dd <- read.table(tree_file, header = T)

  if(branch_num == 0){
    dd$length <- as.numeric(dd$length)
    small_col <- which(dd$length < MIN_BLEN)
    if(length(small_col) > 0){
      dd[small_col, ]$length = 0
    }
    dd <- dd[, c("start", "end", "length")]
  }else if(branch_num == 1){
    dd$nmut <- as.numeric(dd$nmut) * scale_factor
    dd <- dd[, c("start", "end", "nmut")]
  }else{
    stop("branch type not supported!")
  }

  mytree <- make.tree(dd, labels)

  return(mytree)
}


# Get the time differences between patient age at 1st time point for all tips
prepare.tree.age <- function(mytree, time_file){
  # The lengths of tips are at the beginning
  elens = node.depth.edgelength(mytree)

  stime = read.table(time_file, header = F)
  names(stime) = c("sample", "tdiff", "age")

  s1_info = stime[stime$tdiff==0, ][1,]  # samples taken at 1st time point
  age = s1_info$age
  tshift = age - elens[s1_info$sample]  # time difference from born time, LUCA

  # find maximum time difference
  max_tdiff = max(stime$tdiff)

  return(list(tshift = tshift, age = age, max_tdiff = max_tdiff))
}


# Get the tree with bootstrap support value
get.bootstrap.tree <- function(mytree, bstrap_dir, pattern){
  btrees = list()
  cat("Patterns to match bootstrapping trees: ", pattern, "\n")
  files = list.files(path = bstrap_dir, pattern = glob2rx(pattern), recursive = F)
  for(i in 1:length(files)){
    fname = file.path(bstrap_dir, files[i])
    dt = read.table(fname,header = T)
    btree = make.tree(dt, labels)
    btrees[[i]] = btree
  }
  
  # prop.clades calls internally prop.part with the option check.labels = TRUE
  # prop.part counts the number of bipartitions found in a series of trees given as .... If a single tree is passed, the returned object is a list of vectors with the tips descending from each node (i.e., clade compositions indexed by node number).
  clad <- prop.clades(mytree, btrees, rooted = TRUE)
  clad[is.na(clad)] = 0
  mytree$node.label = as.integer(clad * 100 / length(files))
  
  return(mytree)
}


# get the tree with confidence intervals at internal nodes
get.ci.tree <- function(tree_file_nex, bstrap_dir){
  bstrees = list.files(bstrap_dir, pattern = "*.nex")
  ntimes_all = data.frame()
  nbs = length(bstrees)
  for(i in 1:nbs){
    # i = 2
    ftree = file.path(bstrap_dir, bstrees[i])
    btree = read.nexus(ftree)
    
    # get node time
    ndepths = node.depth.edgelength(btree)
    max_depth = max(ndepths)
    ntimes = as.data.frame(max_depth - ndepths)
    colnames(ntimes) = ""
    ntimes_all = rbind(ntimes_all, t(ntimes))
  }
  
  # map Id to node label
  id_labels = c(btree$tip.label, btree$node.label)
  colnames(ntimes_all) = id_labels
  
  # read the tree string
  stree = read_file(tree_file_nex)
  # compute CI interval and append CI to the tree for visualization
  start_inode = btree$Nnode + 2
  root = start_inode
  for(i in start_inode: ncol(ntimes_all)){
    # i = 8
    lbl = id_labels[i]
    #print(lbl)
    times = ntimes_all[, lbl]
    smu = mean(times)
    error = qnorm(0.975) * sd(times) / sqrt(length(times))
    left = smu - error
    right = smu + error
    # [&time_95%_CI={4.50612194676725E-002,1.10259531669597E-001}]
    if(i == start_inode){
      marker = ";"
    }else{
      marker = ":"
    }
    ci_str = paste0(lbl, "[&time_0.95_CI={", left, ",", right, "}]", marker)
    orig_str = paste0(lbl, marker)
    stree = str_replace(stree, orig_str, ci_str)
  }
  
  # write annotated tree to a file
  ftree_ci = str_replace(tree_file, ".txt", ".ci")
  write_lines(stree, ftree_ci)
  
  # read from file and plot CI
  tree_ci = read.mega(ftree_ci)
  lbl_orders = 1:length(labels)
  # ensure labels are ordered by original tree_ci@phylo$tip.label
  tree_ci@phylo$tip.label = labels[match(tree_ci@phylo$tip.label, lbl_orders)]
  
  return(tree_ci)
}


# Add mutation rate for each plot
get.plot.title <- function(mut_rate, dup_rate, del_rate, digit){
  title <- ""
  
  if(mut_rate > 0){
    title <- paste0("mutation rate ", round(mut_rate, digit))
  }
  
  if(dup_rate > 0){
    title <- paste0("duplication rate ", round(dup_rate, digit), "\ndeletion rate ", round(del_rate, digit))
  }
  
  return(title)
}


# annotate tip nodes of a tree 
get.plot.annot <- function(tree, da, p){
  df_labels = data.frame(apeID = 1:nrow(da), Sample = tree$tip.label)
  ds = merge(da, df_labels, by = c("Sample"))
  
  # For IBD data
  if("LocationSpecific" %in% names(da)){
    # color for histology, shape for location
    ds = ds %>% select(node = apeID, HistologySpecific, TP, LocationSpecific)
    ds = ds %>% replace_na(list(LocationSpecific="Unknown", TP = "Unknown", HistologySpecific = "Unknown"))
    p <- p %<+% ds + geom_tippoint(aes(color = HistologySpecific, shape = LocationSpecific), size = 3)
  }
  
  # For BE data
  if("hasWGD" %in% names(da)){
    ds = ds %>% select(node = apeID, hasWGD)
    p <- p %<+% ds + geom_tippoint(aes(shape = as.factor(hasWGD)), size = 3) + theme(legend.position = "none")
  }
  
  # For EAC data
  if("loc" %in% names(da)){
    # color for time, shape for location
    ds = ds %>% select(node = apeID, reltime, loc)
    p <- p %<+% ds + geom_tippoint(aes(color = loc), size = 3) + theme(legend.position = "none")
  }
  
  return(p)
}


plot.tree <- function(tree, title = "", da = data.frame()) {
  p <- ggtree(tree)  #+ geom_rootedge()
  p <- p + geom_tiplab()
  p <- p + geom_text2(aes(subset = !isTip, label = label), hjust = -0.3)
  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge) = c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1) + ggtitle(title)
  
  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }
  
  return(p)
}


# Plot tree with xlim specified to show full tip labels
plot.tree.xlim <- function(tree, title = "", rextra = 20, da = data.frame()) {
  p <- ggtree(tree, size = 0.5, linetype = 1)  #+ geom_rootedge()
  # Add margin to show full name of labels  if (is.na(tree.max))
  tree.max = max(node.depth.edgelength(tree)) + rextra
  p <- p + geom_tiplab(align = TRUE) + theme_tree2() + xlim(NA, tree.max)
  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge) = c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1) + ggtitle(title)
  
  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }
  
  return(p)
}


# Plot tree with xlim specified with age forwards
plot.tree.xlim.age <- function(tree, res_age, title = "", lextra = 3, rextra = 3, da = data.frame()) {
  p <- ggtree(tree)  #, size = 0.5, linetype = 1, + geom_rootedge()

  # Add margin to show full name of labels
  age = res_age$age  # age at 1st sample
   # to display tip text
  tree_max = age + res_age$max_tdiff + rextra
  tree_depth = max(node.depth.edgelength(tree))
  xl = 0
  #xl = res_age$tshift - lextra

  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge) = c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1, nudge_x = res_age$tshift) + ggtitle(title)
  
  # Shift all nodes by the difference between age util the first sample and node time of the first sample
  p$data$x = p$data$x + res_age$tshift
  p <- p + geom_tiplab(align = TRUE) + theme_tree2() + xlim(xl, tree_max) + xlab("Patient age (years)")
  p <- p + geom_text2(aes(subset=!isTip, label=label), hjust=-.3)

  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }

  return(p)
}


plot.tree.bootstrap <- function(tree, title = "", rextra = 20, da = data.frame()){
  p <- ggtree(tree) #+ geom_rootedge()
  # support <- character(length(tree$node.label))
  # #The following three lines define your labeling scheme.
  # support[tree$node.label >= 95] <- "red"
  # support[tree$node.label < 95 & tree$node.label >= 70] <- "pink"
  # support[tree$node.label < 70] <- "blue"
  tree.max = max(node.depth.edgelength(tree)) + rextra
  p <- p + geom_tiplab(align = TRUE) + theme_tree2() + xlim(NA, tree.max)
  p <- p + geom_text2(aes(subset=!isTip, label=label), hjust=-.3, color="red")
  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge)=c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1) + ggtitle(title)
  
  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }
  
  return(p)
}


# Plot bootstrapped tree with x-axis being patient age
# age: used to set limit of x-axis
plot.tree.bootstrap.age <- function(tree, res_age, title = "", lextra = 3, rextra = 3, da = data.frame()){
  p <- ggtree(tree) #+ geom_rootedge()
  # Add margin to show full name of labels

  age = res_age$age  # age at 1st sample
   # to display tip text
  tree_max = age + res_age$max_tdiff + rextra
  tree_depth = max(node.depth.edgelength(tree))
  xl = res_age$tshift - lextra

  # Shift all nodes by the difference between age util the first sample and node time of the first sample
  p$data$x = p$data$x + res_age$tshift
  p <- p + geom_tiplab(align = TRUE) + theme_tree2() + xlim(xl, tree_max) + xlab("Patient age (years)")
  p <- p + geom_text2(aes(subset=!isTip, label=label), hjust=-.3, color="red")
  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge)=c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1, nudge_x = res_age$tshift) + ggtitle(title)

  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }

  return(p)
}


# Plot CI for node time estimation
# read all the trees in the bootstrap folder, bsdir, to compute time CI for original tree
plot.tree.ci.node <- function(tree_ci, time_file, title = "", lextra = 3, rextra = 3, da = data.frame()){
  res_age = prepare.tree.age(tree_ci@phylo, time_file)
  age = res_age$age  # age at 1st sample
  tree_max = age + res_age$max_tdiff + rextra
  tree_depth = max(node.depth.edgelength(tree_ci@phylo))
  xl = res_age$tshift - lextra

  p <- ggtree(tree_ci)
  # Shift all nodes by the difference between age util the first sample and node time of the first sample
  p$data$x = p$data$x + res_age$tshift
  p <- p + geom_tiplab(align = TRUE, offset = 0.5) + theme_tree2() + ggtitle(title)

  # show label of internal nodes
  # p <- p + geom_text2(aes(subset = !isTip, label = label), hjust = -0.3)

  # branch length not properly show after narrow down xaxis
  # edge = data.frame(tree_ci@phylo$edge, edge_num = 1:nrow(tree_ci@phylo$edge), edge_len = tree_ci@phylo$edge.length)
  # colnames(edge)=c("parent", "node", "edge_num", "edge_len")
  # p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1, nudge_x = res_age$tshift)

  p <- p + geom_range('time_0.95_CI', color='red', size=3, alpha=0.5)
  p <- p + xlim(xl, tree_max) + xlab("Patient age (years)")

  if(nrow(da) > 0){
    p = get.plot.annot(tree_ci@phylo, da, p)
  }

  return(p)
}


# Plot a single tree in different format
print.tree <- function(mytree, fout, tree_style, time_file="", tree_file_nex = "", bstrap_dir2 = "", title = "", lextra = 0, rextra = 20, da = data.frame()){

}


option_list = list(
  make_option(c("-f", "--tree_file"), type="character", default="",
              help="dataset file name [default=%default]", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default="",
              help="The name of output file [default=%default]", metavar="character"),
  make_option(c("-a", "--annot_file"), type="character", default="",
              help="The file containing the labels of tip nodes [default=%default]", metavar="character"),
  make_option(c("", "--time_file"), type="character", default="",
              help="The file containing the sampling time information [default=%default]", metavar="character"),
  make_option(c("-d", "--tree_dir"), type="character", default="",
              help="The directory containing all the tree files to plot [default=%default]", metavar="character"),
  make_option(c("-s", "--bstrap_dir"), type="character", default="",
              help="The directory containing all the bootstrapping tree files for tree inference [default=%default]", metavar="character"),
  make_option(c("", "--bstrap_dir2"), type="character", default="",
              help="The directory containing all the bootstrapping tree files when tree topology is fixed [default=%default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="",
              help="The naming pattern of tree files [default=%default]", metavar="character"),
	make_option(c("-b", "--branch_num"), type="integer", default = 0,
              help="The type of values on branches (0: time in year, 1: mutation number) [default=%default]", metavar="integer"),
  make_option(c("", "--scale_factor"), type="numeric", default = 1,
              help="The scale factor to convert branch length unit [default=%default]", metavar="numeric"),
  make_option(c("", "--lextra"), type="numeric", default = 3,
              help="The left margin for the plot [default=%default]", metavar="numeric"),
  make_option(c("", "--rextra"), type="numeric", default = 3,
              help="The left margin for the plot [default=%default]", metavar="numeric"),
  make_option(c("-m", "--mut_rate"), type="numeric", default = 0,
              help="The rate of somatic chromosomal aberrations [default=%default]", metavar="numeric"),
  make_option(c("-u", "--dup_rate"), type="numeric", default = 0,
              help="The segment duplication rate of somatic chromosomal aberrations [default=%default]", metavar="numeric"),
  make_option(c("-e", "--del_rate"), type="numeric", default = 0,
              help="The segment deletion rate of somatic chromosomal aberrations [default=%default]", metavar="numeric"),
  make_option(c("-w", "--with_title"), default = FALSE, action = "store_true", 
              help="Showing title [default=%default]", metavar="logical"),
  make_option(c("", "--title"), type="character", default="",
              help="The title of the plot [default=%default]", metavar="character"),
  # make_option(c("-g", "--use_age"), type="integer", default = 0,
  #             help="Showing x-axis as the real age (0: default (root as beginning time), 1: x-axis as real age of the patient) [default=%default]", metavar="integer"),
  make_option(c("-t", "--plot_type"), type="character", default="single",
              help="The type of plot, including: all (plotting all tree files in a directory), single (plotting a single tree file), bootstrap (plotting a single tree file with bootstrapping support) [default=%default]", metavar="character"),
  make_option(c("-l", "--tree_style"), type="character", default="simple",
              help="The style of tree plot, including: simple (a simple tree with tip labels and branch lengths), xlim (adding xlim to the tree), age (x-axis as real age of the patient), and ci (plotting a single tree file with confidence interval of node ages) [default=%default]", metavar="character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

tree_file = opt$tree_file
out_file = opt$out_file
plot_type = opt$plot_type
tree_dir = opt$tree_dir
pattern = opt$pattern
bstrap_dir = opt$bstrap_dir
bstrap_dir2 = opt$bstrap_dir2
tree_style = opt$tree_style
annot_file = opt$annot_file
branch_num = opt$branch_num
scale_factor = opt$scale_factor
mut_rate = opt$mut_rate
dup_rate = opt$dup_rate
del_rate = opt$del_rate
with_title =  opt$with_title
title = opt$title
time_file = opt$time_file
# cat("Parameters used here:\n")
# cat("tree_file:", tree_file, "\n")
# cat("plot_type:", plot_type, "\n")
# cat("tree_dir:", tree_dir, "\n")
# cat("branch_num:", branch_num, "\n")
# cat("pattern:", pattern, "\n")

if(with_title && title == ""){
  digit = 4
  title = get.plot.title(mut_rate, dup_rate, del_rate, digit)
}
cat("plot title \n", title, "\n")

# add extra lengths along x-axis to show long labels
lextra = opt$lextra
rextra = opt$rextra
if(branch_num == 1){
  rextra = 500
}

if(out_file == ""){
  out_file = get.outfile.name(tree_file)
}

if(file.exists(annot_file)){
  da = read.table(annot_file, header = T, stringsAsFactors = F, sep = "\t")
  if("sample" %in% names(da)){
    da = da %>% dplyr::rename(Sample = sample)
  }else{
    if(!("Sample" %in% names(da))){
      stop("Either \"sample\" or \"Sample\" must be one column of the file!")
    }
  }
  labels = as.character(da$Sample)
}else{
  da = data.frame()
  labels = c()
}


if(plot_type == "all"){
  dir <- tree_dir
  cat(paste0("Plotting all trees in directory ", dir, "\n"))
  if(pattern == ""){
    files <- list.files(dir, "^sim\\-data\\-\\d+\\-tree.txt")
  }else{
    files <- list.files(dir, pattern = glob2rx(pattern))
  }
  #print(files)

  for(f in files){
    cat("running on: ", f, "\n")
    fname = file.path(dir, f)
    tree = get.tree(fname, branch_num, labels, scale_factor)
    print.tree(tree$mytree, tree$fout, tree_style, time_file, tree_file_nex, bstrap_dir2, title, lextraa, rextra)
  }

}else if(plot_type == "single"){
  cat(paste0("Plotting the tree in ", tree_file))

  if(tree_style != "ci"){
    mytree = get.tree(tree_file, branch_num, labels, scale_factor)
    if(tree_style == "simple"){
      p = plot.tree(mytree, title, da)
    }else if(tree_style == "xlim"){
      p = plot.tree.xlim(mytree, title, rextra, da)
    }else if(tree_style == "age"){
      res_age = prepare.tree.age(mytree, time_file)
      p = plot.tree.xlim.age(mytree, res_age, title, lextra, rextra, da)
    }else{
      stop("tree plot style not supported!")
    }
  }else{ # if(tree_style == "ci")
    # cat(paste0("Plotting confidence intervals for the tree in ", tree_file_nex, " with bootstrap trees in folder ", bstrap_dir2, "\n"))
    tree_ci = get.ci.tree(tree_file_nex, bstrap_dir2)
    p = plot.tree.ci.node(tree_ci, out_file, title, lextra, rextra, da)
  }
  
  ggsave(out_file, p, width = 8, height = 5)
  
}else if(plot_type == "bootstrap"){
  cat(paste0("Plotting bootstrap values for the tree in ", tree_file))

  mytree = get.tree(tree_file, out_file, branch_num, labels, scale_factor)
  mytree = get.bootstrap.tree(mytree, bstrap_dir, pattern)

  if(tree_style == "age"){
    res_age = prepare.tree.age(mytree, time_file)
    p = plot.tree.bootstrap.age(mytree, out_file, res_age, title, lextra, rextra, da)
  }else if(tree_style == "ci"){
    # TODO: write bootstrap tree to a nex file 
    tree_ci = get.ci.tree(tree_file_nex, bstrap_dir2)
    p = plot.tree.ci.node(tree_ci, out_file, title, rextra, da)
  }else{
    p = plot.tree.bootstrap(mytree, out_file, title, rextra, da)
  }
  
  ggsave(out_file, p, width = 8, height = 5)
  
}else{
  message("plotting type not supported!")
}
