#!/usr/bin/env Rscript

suppressMessages(library(phangorn))
suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(optparse))

# This script is used to create a set of maximum parsimony trees in NEWICK format,
# serving as initial trees for ML tree building in svtreeml

# convert all CNs into a data frame for tree building
get_phydata <- function(data_cn){
  data_cn %>% unite(chr, seg, col = "pos", sep = "_") %>% spread(pos, cn) %>% select(-sid) -> cns_wide
  # data: For matrices as.phyDat() assumes that the entries each row belongs to one individual (taxa), but for data.frame each column
  data = as.data.frame(t(cns_wide))
  norm = ncol(data)  # normal sample is the last one
  colnames(data) = seq(1, norm)
  dlevel = seq(min(data), max(data))
  phydata <- as.phyDat(data, type='USER', level=dlevel)

  return(phydata)
}


# Convert CN data to breakpoints
get_phydata_bp <- function(data_cn){
  data_cn %>% spread(key = sid, value = cn) %>% select(-c(chr, seg))-> cn_tbl
  bp_list = list()
  for(i in 1:(ncol(cn_tbl))){
    # print(i)
    s1=cn_tbl[,i]
    diff=s1-lag(s1)
    s1_bp=ifelse(diff[-1]==0,0,1)
    bp_list[[i]]=s1_bp
  }

  bp_df = as.data.frame(bp_list)
  names(bp_df) = seq(1,ncol(bp_df))
  dlevel = seq(min(data), max(data))
  phydata <- as.phyDat(data, type='USER', level=dlevel)

  return(phydata)
}


plot_tree <- function(tree){
  p = ggtree(tree)
  p <- p  + geom_tiplab(align = TRUE) + theme_tree2()
  edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge) = c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)
  return(p)
}


build_pratchet <- function(phydata){
  # Build Maximum Parsimony tree
  mptree <- pratchet(phydata)
  mptree <- root(mptree, as.character(ncol(bp_df)), resolve.root=TRUE)
  mptree <- acctran(mptree, phydata)
  plot(mptree)
  # write.nexus(mptree, file=paste0("tree_mp", ".nhs"))
  write.tree(mptree, file=paste0("tree_mp", ".nwk"))
}


pick_top_ntree <- function(uniq_strees, num_select, dir_nwk, format = "NHS"){
  # Sort all trees by score
  scores = c()
  for(i in 1:length(uniq_strees)){
    # print(i)
    stree = uniq_strees[[i]]
    # plot(stree)
    score = attributes(stree)$pscore
    scores = c(scores, score)
  }
  sindex = order(scores)

  # Only take top trees
  num = ifelse(length(uniq_strees) > num_select, num_select, length(uniq_strees))
  sel_trees = list()
  for(i in 1:num){
    j = sindex[i]
    stree = uniq_strees[[j]]
    sel_trees[[i]] = stree
    # score = attributes(stree)$pscore
    if(format == "NWK"){
      write.tree(stree, file=file.path(dir_nwk, paste0("tree_mp_", i, ".nwk")))
    }
  }

  # all.equal.phylo(strees[[1]],strees[[3]])
  # comparePhylo(strees[[1]],strees[[3]])
  if(format == "NHS"){
    write.nexus(sel_trees, file=file.path(dir_nwk, "trees_mp.nhs"), translate = F)
  }
}


# Get input CNs, add normal sample if required
get_cn <- function(file_cn){
  data_cn <- read.table(file_cn)
  names(data_cn) = c("sid", "chr", "seg", "cn")
  # data_cn %>% group_by(sid) %>% count()

  # Check if the last sample is nomal
  nid = length(unique(data_cn$sid))
  data_cn %>% filter(sid == nid) -> sn
  ncn = unique(sn$cn)

  if(length(ncn) == 1 && ncn[1] == 2){
    cat("The last sample is normal\n")
  }else{
    # add a normal sample, use filter() to get segments
    data_cn %>% filter(sid == 1) -> s1
    s1$cn = 2
    s1$sid = length(unique(data_cn$sid)) + 1
    data_cn = rbind(data_cn, s1)
  }

  return(data_cn)
}



option_list = list(
  make_option(c("-f", "--file_cn"), type="character", default="",
              help="input copy number file [default=%default]", metavar="character"),
  make_option(c("-i", "--input_format"), type="integer", default=0,
              help="The input format for tree building (0: copy number; 1: breakpoint)  [default=%default]", metavar="number"),
  make_option(c("-o", "--output_format"), type="integer", default=0,
              help="The output format of selected trees (0: NHS; 1: NWK)  [default=%default]", metavar="number"),
  make_option(c("-d", "--dir_nwk"), type="character", default="",
              help="The directory to store results (NEWICK trees) [default=%default]", metavar="character"),
  make_option(c("-n", "--num_generate"), type="integer", default=100,
              help="The number of trees to generate [default=%default]", metavar="number"),
  make_option(c("-s", "--num_select"), type="integer", default=100,
              help="The maximum number of trees to select [default=%default]", metavar="number")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

file_cn = opt$file_cn
input_format = opt$input_format
num_generate = opt$num_generate
num_select = opt$num_select
dir_nwk = opt$dir_nwk
output_format = ifelse(opt$output_format, "NWK", "NHS")

# file_cn = "D:/data/sveta/test4paper/largetrees/run-5-0.001-0.001-0-1000-1-6/sim-data-100-cn.txt.gz"
# input_format = 0
# num_generate = 100
# num_select = 100
# dir_nwk = "D:/data/sveta/test4paper/largetrees/run-5-0.001-0.001-0-1000-1-6/init_trees"


data_cn = get_cn(file_cn)

if(input_format == 0){
  phydata = get_phydata(data_cn)
}else{
  phydata = get_phydata_bp(data_cn)
}


strees = list()
for(i in seq(1, num_generate)){
  stree = random.addition(phydata, method = "fitch")
  stree <- root(stree, as.character(length(phydata)), resolve.root=TRUE)
  stree <- acctran(stree, phydata)
  # p = plot_tree(stree)
  # print(p)
  strees[[i]] = stree
}

uniq_strees = unique(strees)
cat("There are", length(uniq_strees), "unique trees\n")

dir.create(file.path(dir_nwk), showWarnings = FALSE)
# Remove old files
oldfiles <- dir(path=dir_nwk, pattern="nwk")
if(length(oldfiles) > 0){
  file.remove(file.path(dir_nwk, oldfiles))
}

pick_top_ntree(uniq_strees, num_select, dir_nwk, output_format)
