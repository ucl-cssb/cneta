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
  make_option(c("-n", "--num_tree"), type="integer", default=2000,
              help="The number of trees to generate [default=%default]", metavar="character")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# The file containing copy numbers
file_cn=opt$file_cn
# The directory to store NEWICK trees
dir_nwk=opt$dir_nwk
num_tree=opt$num_tree

dir.create(file.path(dir_nwk), showWarnings = FALSE)
# Remove old files
junk <- dir(path=dir_nwk, pattern="nwk")
# print(junk)
file.remove(file.path(dir_nwk,junk))

data_cn <- read.table(file_cn)
names(data_cn) = c("sid", "chr", "seg", "cn")

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
for(i in seq(1,num_tree))
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


for(i in 1:length(uniq_strees))
{
  print(i)
  stree = uniq_strees[[i]]
  # plot(stree)
  write.tree(stree, file=file.path(dir_nwk, paste0("tree_mp_", i, ".nwk")))
}

# all.equal.phylo(strees[[1]],strees[[3]])
# comparePhylo(strees[[1]],strees[[3]])
