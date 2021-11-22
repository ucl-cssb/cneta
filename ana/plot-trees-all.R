#!/usr/bin/env Rscript

suppressMessages(library(ggtree))
suppressMessages(library(treeio))
suppressMessages(library(ape))
suppressMessages(library(tidyverse))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(optparse))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
# suppressMessages(library(aplot))

# This script is used to plot phylogenetic trees:
# 1) plotting a single tree file with tips as 1 to n or with provided labels
# 2) plotting a single tree file with bootstrapping support
# 3) plotting all tree files in a directory


theme1 = theme(legend.position = "none",
               #strip.text.x = element_blank(),
               #strip.text.y = element_blank(),
               strip.text.y.left = element_text(size=6, angle = 0),
               strip.background = element_blank(),
               axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.ticks.y = element_blank(),
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank(),
               panel.spacing.y = unit(0, "lines"),
               panel.spacing.x = unit(0, "lines"),
               panel.border = element_rect(color = "#d8d8d8", fill = NA, size = .2),
               panel.background = element_rect(fill = "#f0f0f0")
)

# cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0","#FCAE91", "#B9574E", "#76000D", "#8B0000", "#000000")
# Max CN to show in heatmap
MAX_CN = 8
# For absolute CN 0 to 8
cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0", '#fdd49e', '#fdbb84','#fc8d59','#ef6548','#d7301f','#990000')
# For relative CN -4 to 4
cn_colors2 = c('#08519c','#3182bd', '#6baed6', '#9ecae1', "#f0f0f0", '#fdd49e','#fc8d59','#d7301f','#990000')


MIN_BLEN = 1e-3
TIP_OFFSET = 0.5
BLEN_DIGIT = 2

# d is a data frame with 3 columns: start, end, branch length
# Leaves must be encoded from 1 to nleaf
make.tree <- function(d, labels = c()) {
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
  mytree$edge.length <- d[, 3]
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

  # cat("\n\nrunning over", stub, fout, "\n", sep = "\t")

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


get_cyto_band <- function(fcyto){
  cyto = read.table(fcyto, header=FALSE)
  names(cyto) <- c("chrom", "start", "end", "cytoband", "stain")
  
  # Find the boundary between two arms
  arm_boudaries = data.frame()
  chr_ends = data.frame()
  for (i in 1:(nrow(cyto)-1)) {
    chr1 = cyto[i,]$chrom
    arm1 = substr(cyto[i,]$cytoband, 1, 1)
    chr2 = cyto[i+1,]$chrom
    arm2 = substr(cyto[i+1,]$cytoband, 1, 1)
    
    if(chr1 == chr2 && arm1 != arm2){
      item = data.frame(chrom = chr1, boundary = cyto[i+1,]$start)
      arm_boudaries = rbind(arm_boudaries, item)
      # print(paste("chrom", chr1, "arm boundary", cyto[i+1,]$start))
    }
    if(chr1 != chr2){
      item = data.frame(chrom=chr1, cend = cyto[i,]$end)
      chr_ends = rbind(chr_ends, item)
    }
    if(i == nrow(cyto) - 1){
      item = data.frame(chrom=chr2, cend = cyto[i+1,]$end)
      chr_ends = rbind(chr_ends, item)      
    }
  }
  
  merge(chr_ends, arm_boudaries, by="chrom") -> chr_end_arm
  
  return(chr_end_arm)
}


get_hg_ends <- function(fcyto){
  chr_end_arm_hg = get_cyto_band(fcyto)
  
  chr_end_arm_hg %>% filter(chrom!="chrX" & chrom!="chrY") -> chr_end_arm_hg_short
  # Use cumsum to convert coordinates in the data
  chr_end_arm_hg_short$chrID = gsub("chr", "", chr_end_arm_hg_short$chrom)
  chr_end_arm_hg_short %>% arrange(as.numeric(chrID)) -> chr_end_arm_hg_short
  csize = cumsum(as.numeric(chr_end_arm_hg_short$cend)+1)
  chr_end_arm_hg_short$csize = c(0, csize[1:length(csize)-1])
  chr_end_arm_hg_short %>% select(chromosome = chrID, cend, csize, boundary) -> chr_end_arm_hg_short
  
  return(chr_end_arm_hg_short)
}


# add missing segments 
get_normal_segs <- function(d_withpos, chr_end_arm){
  samples = data.frame(sample = unique(d_withpos$sample))
  
  d_noindex = d_withpos  %>% select(-index)
  
  # get the end of intermediate intervals first
  df_invar = d_noindex %>% dplyr::mutate(end = lead(start) - 1) 
  df_invar$start = d_withpos$end + 1
  df_invar = df_invar %>% filter(end > start)
  df_invar$cn = 2
  
  # Update interval at the end of each chr
  # Find the current end in the data for each chr
  d_withpos %>% group_by(chromosome) %>% dplyr::summarise(start = max(end)) -> pos_max
  seg_max = merge(pos_max, chr_end_arm, by = c("chromosome")) %>% dplyr::mutate(start = start + 1) %>% select(chromosome, start, end = cend) %>% filter(start < end)
  cnT_max = merge(seg_max, samples, all = T) %>% arrange(sample, chromosome, start)
  cnT_max$cn = 2
  
  
  # add intervals at the beginning of each chr
  d_withpos %>% group_by(chromosome) %>% dplyr::summarise(end = min(start)) -> pos_min
  seg_min = pos_min %>% dplyr::mutate(start = 1, end = end - 1) %>% filter(start < end)
  cnT_min = merge(seg_min, samples, all = T) %>% arrange(sample, chromosome, start)
  cnT_min$cn = 2 
  
  # Remove interval at the end which has wrong boundaries
  max_ends = pos_max$start + 1
  min_ends = pos_min$end - 1
  row2del = df_invar %>% filter(start %in% max_ends | end %in% min_ends)
  df_invar_mid = anti_join(df_invar, row2del)
  
  # add missing full chr
  chrs = unique(d_withpos$chromosome)
  all_chrs = unique(chr_end_arm$chromosome)
  miss_chrs = setdiff(all_chrs, chrs)
  seg_chm = data.frame()
  for(ch in miss_chrs){
    seg_df = chr_end_arm %>% filter(chromosome == ch) %>% dplyr::mutate(chromosome = as.numeric(chromosome), start = 1) %>% select(chromosome, start, end = cend) 
    cnT_df = merge(seg_df, samples, all = T) 
    cnT_df$cn = 2 
    seg_chm = rbind(seg_chm, cnT_df)
  }
  
  # Add Seg for all regions
  cn_all = rbind(d_noindex, cnT_min, df_invar_mid, cnT_max, seg_chm) %>% arrange(sample, chromosome, start) %>% group_by(sample) %>% dplyr::mutate(index = row_number()) %>% select(sample, chromosome, index, start, end, cn)
  
  return(cn_all)
}


# combine CNP with position information to get a file with the format: sample, chrom, start, end, cn
# use chr and site index to bind the two datasets
get.cn.data.by.pos <- function(cn_file, pos_file, seg_file, cyto_file, labels, ordered_nodes, add_normal = T){
  # reconstructed states
  dans <- read.table(cn_file, header=FALSE)
  names(dans) <- c("sample", "chromosome", "index", "cn")
  ans_nodes = unique(dans$sample) %>% unlist()
  ans_labels = data.frame(sample=ans_nodes, name = ans_nodes)
  
  # original input, used to get sample names and original positions
  dpos <- read.table(pos_file, header=FALSE)
  names(dpos) <- c("sample", "chromosome", "index", "cn", "start", "end")
  
  samples = unique(dpos$sample) %>% unlist()
  tip_labels = data.frame(sample = 1:length(labels), name = labels)
  
  # extract positions for all sites
  dpos %>% select(chromosome, index, start, end) %>% unique() -> site_pos
  
  # read segment file for original input as invariable bins are excluded and consecutive bins may be merged
  if(seg_file == ""){
    stop("a segment file must be provided!")
  }
  
  dseg <- read.table(seg_file, header=FALSE)
  segs = dseg[, c(1,4,5)]
  names(segs) = c("chromosome", "index", "index2")
  # first get positions for segment start
  segs_pos = merge(segs, site_pos, by = c("chromosome", "index"), all.x = T)
  # then get positions for segment end
  segs_pos2 = merge(segs_pos, site_pos, by.x = c("chromosome", "index2"), by = c("chromosome", "index"), all.x = T)
  segs_pos_all = segs_pos2 %>% select(chromosome, start = start.x, end = end.y) %>% arrange(start, end) %>% group_by(chromosome)  %>%  mutate(index = 1:n())
  
  # convert segment into cn matrix to have consistent segments in both input and reconstructed states
  # wide format to long format 
  dseg_cn = dseg[, -c(3,4,5)]
  # add "s" to avoid taking first n columns when using gather
  names(dseg_cn) = c("chromosome", "index", paste0("s", all_of(samples)))
  dseg_cn_long = dseg_cn %>% gather(paste0("s", all_of(samples)), key = "sample", value = "cn") %>% mutate(sample = str_replace(sample, "s","")) %>% group_by(sample, chromosome)  %>%  mutate(index = 1:n()) %>% select(sample, chromosome, index, cn)
  dseg_cn_long$sample = as.numeric(dseg_cn_long$sample)
  
  # add normal samples
  if(add_normal){
    # add a normal sample, use filter() to get segments
    dseg_cn_long %>% filter(sample == 1) -> s1
    s1$cn = 2
    s1$sample = length(unique(dseg_cn_long$sample)) + 1
    dseg_cn_long = rbind(dseg_cn_long, s1)
  }

  d_binded = rbind(dans, dseg_cn_long)
  
  # get reconstructed states with position
  d_all = merge(d_binded, segs_pos_all, by = c("chromosome", "index"), all.x = T) %>% select(all_of(names(dpos))) %>% arrange(sample) 
  # }else{ # using initial bins
  #   # add normal samples
  #   if(add_normal){
  #     # add a normal sample, use filter() to get segments
  #     dpos %>% filter(sample == 1) -> s1
  #     s1$cn = 2
  #     s1$sample = length(unique(dpos$sample)) + 1
  #     dpos = rbind(dpos, s1)
  #   }
  # 
  #   d_ans = merge(dans, site_pos, by = c("chromosome", "index"), all.x = T) %>% select(all_of(names(dpos))) %>% arrange(sample)
  #   d_all = rbind(dpos_rname, d_ans)
  # }
  # 
  # add missing normal segments
  chr_end_arm = get_hg_ends(cyto_file)
  d_all = get_normal_segs(d_all, chr_end_arm) 
  
  # replace node IDs with sample names 
  node_labels = rbind(ans_labels, tip_labels)
  d_all = merge(d_all, node_labels, by = c("sample"), all.x = T) %>% select(-sample) %>% select(sample = name, chromosome, index, start, end, cn)
  d_all$sample = factor(d_all$sample, levels = ordered_nodes)
  
  # d_all %>% group_by(chromosome , sample) %>% tally() %>% select(-sample) %>% unique() %>% print()
  
  d_seg = d_all %>% select(sample, chrom = chromosome, start, end, cn)
  
  return(d_seg)
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
get.ci.tree <- function(tree_file_nex, bstrap_dir, labels, has_bstrap = F){
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
  # all the bootstrapping trees have the same topology here
  id_labels = c(btree$tip.label, btree$node.label)
  colnames(ntimes_all) = id_labels

  # read the tree string
  stree = read_file(tree_file_nex)
  ape_tree = read.nexus(tree_file_nex)

  # compute CI interval and append CI to the tree for visualization
  start_inode = btree$Nnode + 2
  root = start_inode

  for(i in start_inode: ncol(ntimes_all)){
    # i = 8
    lbl = id_labels[i]
    # print(i)
    # print(lbl)
    times = ntimes_all[, lbl]
    
    # smu = mean(times)
    # error = qnorm(0.975) * sd(times) / sqrt(length(times))
    # left = smu - error
    # right = smu + error
    sorted_times = sort(times)
    ci95 = quantile(sorted_times, probs = c(0.025, 0.975))
    left = ci95[1]
    right = ci95[2]
    
    # [&time_95%_CI={4.50612194676725E-002,1.10259531669597E-001}]
    if(i == start_inode){
      marker = ";"
      blen = ""
    }else{
      marker = ":"
      # use branch length to avoid mismatching
      # find the edge ID with the node
      eid = which(ape_tree$edge[, 2]==i)
      blen = ape_tree$edge.length[eid]
    }

    if(has_bstrap){
      # when the tree has bootstrap support value, the node label is the support value
      # assume: all internal branches have different lengths
      # find the index of current node
      nid = which(id_labels == lbl)
      bsval = ape_tree$node.label[nid - ape_tree$Nnode - 1]
      orig_str = paste0(bsval, marker, blen)
      ci_str = paste0(lbl, "[&time_0.95_CI={", left, ",", right, "},bootstrap=", bsval, "]", marker, blen)
    }else{
      orig_str = paste0(lbl, marker, blen)
      ci_str = paste0(lbl, "[&time_0.95_CI={", left, ",", right, "}]", marker, blen)
    }
    # print(orig_str)
    # print(ci_str)
    stree = str_replace(stree, orig_str, ci_str)
  }

  # write annotated tree to a file
  ftree_ci = str_replace(tree_file, ".txt", ".ci.nex")
  #print(stree)
  write_lines(stree, ftree_ci)

  # read from file and plot CI
  tree_ci = read.mega(ftree_ci)

  if((!has_bstrap) && length(labels) > 0){
    lbl_orders = 1:length(labels)
    # ensure labels are ordered by original tree_ci@phylo$tip.label
    tree_ci@phylo$tip.label = labels[match(tree_ci@phylo$tip.label, lbl_orders)]
  }
  
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

  #Â For BE data
  if("hasWGD" %in% names(da)){
    ds = ds %>% select(node = apeID, hasWGD)
    p <- p %<+% ds + geom_tippoint(aes(shape = as.factor(hasWGD)), size = 3) + theme(legend.position = "none")
  }

  # For EAC data
  if("loc" %in% names(da)){
    # color for time, shape for location
    ds = ds %>% select(node = apeID, reltime, loc)
    p <- p %<+% ds + geom_tippoint(aes(color = as.factor(loc)), size = 3) + theme(legend.position = "none")
  }

  return(p)
}


plot.tree <- function(tree, title = "", da = data.frame()) {
  p <- ggtree(tree)  #+ geom_rootedge()
  p <- p + geom_tiplab()
  p <- p + geom_text2(aes(subset = !isTip, label = label), hjust = -0.3)

  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge) = c("parent", "node", "edge_num", "edge_len")
  edge$edge_len = round(edge$edge_len, BLEN_DIGIT)

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
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + theme_tree2() + xlim(NA, tree.max)

  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge) = c("parent", "node", "edge_num", "edge_len")
  edge$edge_len = round(edge$edge_len, BLEN_DIGIT)

  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1) + ggtitle(title)

  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }

  return(p)
}


# Plot tree with xlim specified with age forwards
plot.tree.xlim.age <- function(tree, time_file, title = "", lextra = 3, rextra = 3, da = data.frame()) {
  res_age = prepare.tree.age(tree, time_file)

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
  edge$edge_len = round(edge$edge_len, BLEN_DIGIT)

  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1, nudge_x = res_age$tshift) + ggtitle(title)

  # Shift all nodes by the difference between age util the first sample and node time of the first sample
  p$data$x = p$data$x + res_age$tshift
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + theme_tree2() + xlim(xl, tree_max) + xlab("Patient age (years)")
  p <- p + geom_text2(aes(subset=!isTip, label=label), hjust=-.3)

  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }

  return(p)
}


plot.tree.bootstrap <- function(tree, title = "", rextra = 20, da = data.frame()){
  p <- ggtree(tree)
  #+ geom_rootedge()

  # support <- character(length(tree$node.label))
  # #The following three lines define your labeling scheme.
  # support[tree$node.label >= 95] <- "red"
  # support[tree$node.label < 95 & tree$node.label >= 70] <- "pink"
  # support[tree$node.label < 70] <- "blue"

  tree.max = max(node.depth.edgelength(tree)) + rextra
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + theme_tree2() + xlim(NA, tree.max)
  p <- p + geom_text2(aes(subset=!isTip, label=label), hjust=-.3, color="red")

  # edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  # colnames(edge)=c("parent", "node", "edge_num", "edge_len")
  # p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)

  p <- p + ggtitle(title) + xlab("Size of copy number alteration (Mbp)")

  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }

  return(p)
}


# Plot bootstrapped tree with x-axis being patient age
# age: used to set limit of x-axis
plot.tree.bootstrap.age <- function(tree, time_file, title = "", lextra = 3, rextra = 3, da = data.frame()){
  res_age = prepare.tree.age(tree, time_file)

  p <- ggtree(tree)

  # Add margin to show full name of labels
  age = res_age$age  # age at 1st sample
   # to display tip text
  tree_max = age + res_age$max_tdiff + rextra
  tree_depth = max(node.depth.edgelength(tree))
  xl = res_age$tshift - lextra

  # Shift all nodes by the difference between age util the first sample and node time of the first sample
  p$data$x = p$data$x + res_age$tshift
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + theme_tree2() + xlim(xl, tree_max) + xlab("Patient age (years)")
  p <- p + geom_text2(aes(subset=!isTip, label=label), hjust=-.3, color="red")

  # edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  # colnames(edge)=c("parent", "node", "edge_num", "edge_len")
  # p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1, nudge_x = res_age$tshift)

  p <- p + ggtitle(title)

  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }

  return(p)
}


# Plot CI for node time estimation
# read all the trees in the bootstrap folder, bsdir, to compute time CI for original tree
plot.tree.ci.node <- function(tree_ci, time_file, title = "", lextra = 3, rextra = 3, da = data.frame(), has_bstrap = F, has_inode_label = F){
  res_age = prepare.tree.age(tree_ci@phylo, time_file)
  age = res_age$age  # age at 1st sample
  tree_max = age + res_age$max_tdiff + rextra
  tree_depth = max(node.depth.edgelength(tree_ci@phylo))
  xl = res_age$tshift - lextra

  p <- ggtree(tree_ci)
  # Shift all nodes by the difference between age util the first sample and node time of the first sample
  p$data$x = p$data$x + res_age$tshift
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + theme_tree2() + ggtitle(title)
  
  if(has_inode_label){
    p <- p + geom_text2(aes(subset = !isTip, label = label), hjust = -0.3)
  }
 
  # show label of internal nodes
  if(has_bstrap){
    p <- p + geom_label_repel(aes(label = bootstrap, fill = bootstrap)) + scale_fill_viridis_c(alpha = 0.5)
  }

  # branch length not properly show after narrow down xaxis
  # edge = data.frame(tree_ci@phylo$edge, edge_num = 1:nrow(tree_ci@phylo$edge), edge_len = tree_ci@phylo$edge.length)
  # colnames(edge)=c("parent", "node", "edge_num", "edge_len")
  # p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1, nudge_x = res_age$tshift)

  p <- p + geom_range('time_0.95_CI', color='red', size=3, alpha = 0.3)
  p <- p + xlim(xl, tree_max) + xlab("Patient age (years)")

  if(nrow(da) > 0){
     p = get.plot.annot(tree_ci@phylo, da, p)
  }

  return(p)
}


get.cn.matrix <- function(cn_file){
  data_cn = read.table(cn_file)
  names(data_cn) = c("sid", "chr", "seg", "cn")
  data_cn %>% unite(chr, seg, col = "pos", sep = "_") %>% spread(pos, cn) %>% select(-sid) -> cns_wide
  
  return(cns_wide)
}


# Plot CNPs of multiple samples as a heatmap for easy comparison
# Format of d_seg: sample, chrom, start, end, cn, ploidy
plot.cn.heatmap <- function(d_seg, main, type="absolute", theme = theme1, cn_colors = cn_colors1, allele_specific = F){
  d_seg$cn = round(d_seg$cn)
  d_seg$chrom = factor(d_seg$chrom, levels=paste("",c(c(1:22)),sep=""))
  d_seg$pos = (d_seg$end + d_seg$start) / 2
  d_seg$width = (d_seg$pos - d_seg$start) * 2
  
  if(type=="absolute"){
    print("Plot absolute copy number")
    d_seg %>% dplyr::mutate(cn=if_else(cn > MAX_CN, MAX_CN, cn)) -> d_seg
    cn_vals = c("0", "1", "2", "3", "4", "5", "6", "7", "8")
  }else{
    print("Plot relative copy number")
    d_seg %>% dplyr::mutate(cn=if_else(cn > max_rcn, max_rcn, cn)) %>% dplyr::mutate(cn=if_else(cn < min_rcn, min_rcn, cn)) -> d_seg
    cn_vals = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4")
  }
  
  d_seg$cn = as.factor(d_seg$cn)
  # unique(d_seg$cn)
  
  p <- ggplot(data=d_seg, aes(x=pos, y=sample, width=width)) +
    geom_tile(aes(fill=cn)) +
    facet_grid(sample ~ chrom, scales="free", space = "free", switch='both') +
    scale_color_manual(values = cn_colors, limits=cn_vals) +
    scale_fill_manual(values = cn_colors, limits=cn_vals) +
    scale_y_discrete(breaks = unique(d_seg$sample), expand = c(0,0))+
    scale_x_discrete(expand = c(0,0)) + theme
  p = p + ggtitle(main)
  
  return(p)
}


# Plot a single tree in different format
print.single.tree <- function(mytree, tree_style, time_file="", title = "", lextra = 0, rextra = 20, da = data.frame()){
  if(tree_style == "simple"){
    p = plot.tree(mytree, title, da)
  }else if(tree_style == "xlim"){
    p = plot.tree.xlim(mytree, title, rextra, da)
  }else if(tree_style == "age"){
    if(time_file == ""){
      stop("The file containing the sampling time information is not provided!")
    }
    p = plot.tree.xlim.age(mytree, time_file, title, lextra, rextra, da)
  }else{
    stop("tree plot style not supported!")
  }
  return(p)
}


option_list = list(
  make_option(c("-f", "--tree_file"), type="character", default="",
              help="tree file name in TSV format [default=%default]", metavar="character"),
  make_option(c("", "--tree_file_nex"), type="character", default="",
              help="tree file name in NEXUS format [default=%default]", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default="",
              help="The name of output file [default=%default]", metavar="character"),
  make_option(c("-a", "--annot_file"), type="character", default="",
              help="The file containing the labels of tip nodes [default=%default]", metavar="character"),
  make_option(c("", "--time_file"), type="character", default="",
              help="The file containing the sampling time information [default=%default]", metavar="character"),
  make_option(c("", "--cn_file"), type="character", default="",
              help="The file which contains the copy numbers of internal nodes [default=%default]", metavar="character"),
  make_option(c("", "--pos_file"), type="character", default="",
              help="The file which contains the positions of each site along with the copy numbers of tip nodes [default=%default]", metavar="character"),
  make_option(c("", "--seg_file"), type="character", default="",
              help="The file which contains the segment start and end ID in the input copy number file [default=%default]", metavar="character"),
  make_option(c("", "--cyto_file"), type="character", default="",
              help="The file which contains the chromosome boundaries in human reference genome (e.g. hg19) [default=%default]", metavar="character"),
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
  make_option(c("", "--with_cn"), default = FALSE, action = "store_true",
              help="Plotting copy numbers [default=%default]", metavar="logical"),
  make_option(c("", "--title"), type="character", default="",
              help="The title of the plot [default=%default]", metavar="character"),
  make_option(c("-t", "--plot_type"), type="character", default="single",
              help="The type of plot, including: all (plotting all tree files in a directory), single (plotting a single tree file), bootstrap (plotting a single tree file with bootstrapping support) [default=%default]", metavar="character"),
  make_option(c("-l", "--tree_style"), type="character", default="simple",
              help="The style of tree plot, including: simple (a simple tree with tip labels and branch lengths), xlim (adding xlim to the tree), age (x-axis as real age of the patient), and ci (plotting a single tree file with confidence interval of node ages) [default=%default]", metavar="character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

tree_file = opt$tree_file
tree_file_nex = opt$tree_file_nex
out_file = opt$out_file
time_file = opt$time_file
cn_file = opt$cn_file
pos_file = opt$pos_file
seg_file = opt$seg_file
cyto_file = opt$cyto_file
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
with_title = opt$with_title
with_cn = opt$with_cn
title = opt$title

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

  for(f in files){
    fname = file.path(dir, f)
    cat("running on: ", fname, "\n")
    out_file = get.outfile.name(fname)
    mytree = get.tree(fname, branch_num, labels, scale_factor)
    p = print.single.tree(mytree, tree_style, time_file, title, lextra, rextra, da)
    ggsave(out_file, p, width = 8, height = 5)
  }

}else if(plot_type == "single"){
  cat(paste0("Plotting the tree in ", tree_file, "\n"))
  if(out_file == ""){
    out_file = get.outfile.name(tree_file)
  }

  if(tree_style != "ci"){
    mytree = get.tree(tree_file, branch_num, labels, scale_factor)
    p = print.single.tree(mytree, tree_style, time_file, title, lextra, rextra, da)
  }else{
    # cat(paste0("Plotting confidence intervals for the tree in ", tree_file_nex, " with bootstrap trees in folder ", bstrap_dir2, "\n"))
    tree_ci = get.ci.tree(tree_file_nex, bstrap_dir2, labels)
    p = plot.tree.ci.node(tree_ci, time_file, title, lextra, rextra, da)
  }

  ggsave(out_file, p, width = 8, height = 5)

}else if(plot_type == "bootstrap"){
  cat(paste0("Plotting bootstrap values for the tree in ", tree_file, "\n"))

  # when the tree is read from txt file, the internal node labels follow the same order as input CNs
  mytree = get.tree(tree_file, branch_num, labels, scale_factor)
  mytree = get.bootstrap.tree(mytree, bstrap_dir, pattern)

  if(tree_style == "age"){
    if(time_file == ""){
      stop("The file containing the sampling time information is not provided!")
    }
    p = plot.tree.bootstrap.age(mytree, time_file, title, lextra, rextra, da)
  }else if(tree_style == "ci"){
    if(time_file == ""){
      stop("The file containing the sampling time information is not provided!")
    }
    # write bootstrap tree to a nex file
    fbs = str_replace(tree_file, ".txt", ".bstrap.nex")
    write.nexus(mytree, file = fbs)
    # the tree has correct tip labels after using write.nexus
    tree_ci = get.ci.tree(fbs, bstrap_dir2, labels, T)
    
    p = plot.tree.ci.node(tree_ci, time_file, title, lextra, rextra, da, T, T)
    
    if(with_cn){
      d = fortify(tree_ci@phylo)
      ordered_nodes = d$label[order(d$y, decreasing = T)]
      d_seg = get.cn.data.by.pos(cn_file, pos_file, seg_file, cyto_file, labels, ordered_nodes, T)
      # get the node order of the tree and reorder heatmap
      phmap = plot.cn.heatmap(d_seg, "")
      
      pc = ggarrange(p, phmap, nrow = 1, widths = c(6.5, 9.5))  
      p#pc = phmap %>% insert_left(p, width = 0.6)
      ggsave(out_file, pc, width = 16, height = 5)
    }else{
      ggsave(out_file, p, width = 8, height = 5)
    }
    
  }else{
    p = plot.tree.bootstrap(mytree, title, rextra, da)
    ggsave(out_file, p, width = 8, height = 5)
  }

 

}else{
  message("plotting type not supported!")
}
