library(phangorn)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(optparse)

# This script is used to create a set of maximum parsimony trees in NEWICK format

option_list = list(
  make_option(c("-f", "--file_cn"), type="character", default="",
              help="input copy number file [default=%default]", metavar="character"),
  make_option(c("-d", "--dir_nwk"), type="character", default="",
              help="The directory to store results [default=%default]", metavar="character"),
  make_option(c("-n", "--num_generate"), type="integer", default=2000,
              help="The number of trees to generate [default=%default]", metavar="number"),
  make_option(c("-s", "--num_select"), type="integer", default=100,
              help="The maximum number of trees to select [default=%default]", metavar="number")  
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# The file containing copy numbers
file_cn=opt$file_cn
# The directory to store NEWICK trees
dir_nwk=opt$dir_nwk
num_generate=opt$num_generate
num_select=opt$num_select


# file_cn="D:/Gdrive/git/sveta/example/391forSveta-cn.txt"
# dir_nwk="D:/Gdrive/git/sveta/example/itrees_391_nwk"
# num_generate=100
# num_select=10

dir.create(file.path(dir_nwk), showWarnings = FALSE)
# Remove old files
junk <- dir(path=dir_nwk, pattern="nwk")
# print(junk)
file.remove(file.path(dir_nwk,junk))

data_cn <- read.table(file_cn)
names(data_cn) = c("sid", "chr", "seg", "cn")
# data_cn %>% group_by(sid) %>% count()
# Check if the last sample is nomal
nid=length(unique(data_cn$sid)) 
data_cn %>% filter(sid==nid) -> sn
ncn=unique(sn$cn)

if(length(ncn)==1 && ncn[1]==2){
  cat("The last sample is normal")
}else{
  data_cn %>% filter(sid==1) -> s1
  s1$cn=2
  s1$sid=length(unique(data_cn$sid))+1
  data_cn=rbind(data_cn, s1)
}

# nid=length(unique(data_cn$sid)) 
# data_cn %>% filter(sid==nid) -> sn
# unique(sn$cn)

# Convert data to breakpoints
data_cn %>% spread(key = sid, value = cn) %>% select(-c(chr, seg))-> cn_tbl
bp_list=list()
for(i in 1:(ncol(cn_tbl))){
  # print(i)
  s1=cn_tbl[,i]
  diff=s1-lag(s1)
  s1_bp=ifelse(diff[-1]==0,0,1)
  bp_list[[i]]=s1_bp
}

bp_df = as.data.frame(bp_list)
names(bp_df)=seq(1,ncol(bp_df))
# unique(bp_df)
phydata_cnbp <- as.phyDat(bp_df, type='USER', level=(c(0, 1)))


# Build Maximum Parsimony tree
# MPdata_cnbp <- pratchet(phydata_cnbp, start=NULL, k=20, trace=1, all=TRUE, method="fitch")
# MPdata_cnbp <- root(MPdata_cnbp, as.character(ncol(bp_df)), resolve.root=TRUE)
# MPdata_cnbp <- acctran(MPdata_cnbp, phydata_cnbp)
# plot(MPdata_cnbp)
# # write.nexus(MPdata_cnbp, file=paste0("tree_mp", ".nhs"))
# write.tree(MPdata_cnbp, file=paste0("tree_mp", ".nwk"))
#
# tree=MPdata_cnbp
# p = ggtree(tree)
# p <- p  + geom_tiplab(align = TRUE) + theme_tree2()
# edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
# colnames(edge) = c("parent", "node", "edge_num", "edge_len")
# p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)
# print(p)


strees = list()
for(i in seq(1,num_generate))
{
  # set.seed(i)
  stree=random.addition(phydata_cnbp, method = "fitch")
  stree <- root(stree, as.character(ncol(bp_df)), resolve.root=TRUE)
  stree <- acctran(stree, phydata_cnbp)
  # plot(stree)
  strees[[i]] = stree
  # write.tree(stree, file=file.path(dir, paste0("tree_mp_", i, ".nwk")))
}
uniq_strees = unique.multiPhylo(strees, use.edge.length = TRUE)
# length(uniq_strees)

# Sort all trees by score
scores=c()
for(i in 1:length(uniq_strees))
{
  # print(i)
  stree = uniq_strees[[i]]
  score = attributes(stree)$pscore
  scores = c(scores, score)
  # print(score)
  # plot(stree)
  #write.tree(stree, file=file.path(dir_nwk, paste0("tree_mp_", i, ".nwk")))
}
sindex = order(scores)

# Only take top trees
num = ifelse(length(uniq_strees)>num_select, num_select, length(uniq_strees))
for(i in 1:num)
{
  j = sindex[i]
  # print(j)
  stree = uniq_strees[[j]]
  score = attributes(stree)$pscore
  # print(score)
  # plot(stree)
  write.tree(stree, file=file.path(dir_nwk, paste0("tree_mp_", i, ".nwk")))
}

# all.equal.phylo(strees[[1]],strees[[3]])
# comparePhylo(strees[[1]],strees[[3]])
