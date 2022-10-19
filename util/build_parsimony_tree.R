#!/usr/bin/env Rscript

suppressMessages(library(phangorn))
suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(optparse))

# This script is used to create a set of maximum parsimony trees in NEWICK format,
# serving as initial trees for ML tree building in cnetml

# convert all CNs into a data frame for tree building
get_phydata <- function(cns_wide) {
  # data: For matrices as.phyDat() assumes that the entries each row belongs to one individual (taxa), but for data.frame each column
  data <- as.data.frame(t(cns_wide))
  norm <- ncol(data) # normal sample is the last one
  colnames(data) <- seq(1, norm)
  dlevel <- seq(min(data), max(data))
  phydata <- as.phyDat(data, type = "USER", level = dlevel)

  return(phydata)
}


# Convert CN data to breakpoints
get_phydata_bp <- function(data_cn) {
  data_cn %>%
    spread(key = sid, value = cn) %>%
    select(-c(chr, seg)) -> cn_tbl
  bp_list <- list()
  for (i in 1:(ncol(cn_tbl))) {
    # print(i)
    s1 <- cn_tbl[, i]
    diff <- s1 - lag(s1)
    s1_bp <- ifelse(diff[-1] == 0, 0, 1)
    bp_list[[i]] <- s1_bp
  }

  bp_df <- as.data.frame(bp_list)
  names(bp_df) <- seq(1, ncol(bp_df))
  dlevel <- seq(min(data), max(data))
  phydata <- as.phyDat(data, type = "USER", level = dlevel)

  return(phydata)
}


get_bootstrap_phydata <- function(data_cn, is_haplotype_specific) {
  # the column order is changed to be alphabetical
  data_cn %>%
    unite(chr, seg, col = "pos", sep = "_") %>%
    spread(pos, cn) %>%
    select(-sid) -> cns_wide
  # resample columns to get the same dimension

  # when the data is haplotype specific, need to ensure both haplotypes of a site are chosen
  if(is_haplotype_specific){
    nsite <- ncol(cns_wide) / 2
    sites = names(cns_wide)
    uniq_sites = sites[!grepl("-", sites)]
    nsite = length(uniq_sites)
    bs_idx <- sample(1:nsite, nsite, replace = T)
    sel_sites = uniq_sites[bs_idx]
    sites2 = paste0(sel_sites, "-2")
    cns_bs <- cns_wide[, c(sel_sites, sites2)]   
  }else{
    nsite <- ncol(cns_wide)
    bs_idx <- sample(1:nsite, nsite, replace = T)
    cns_bs <- cns_wide[, bs_idx]   
  }

  # build step addition trees for bootstrap sample
  phydata <- get_phydata(cns_bs)

  # convert data into format required by cnets  c("sid", "chr", "seg", "cn")
  # check the convertion is right by converting cns_wide back to original data
  # cns_wide$sid = 1:nrow(cns_wide)
  # pos_names = setdiff(names(cns_wide), c("sid"))
  # cns_wide %>% gather(all_of(pos_names), key = "pos", value = "cn") %>% separate(pos, sep="_", into = c("chr", "seg")) %>% mutate(chr = as.integer(chr), seg = as.integer(seg))%>% arrange(sid, chr, seg) -> dcn

  return(list(phydata = phydata, cns_bs = cns_bs))
}


# output the copy numbers corresponding to bootstrap samples
get_bootstrap_cn <- function(cns_bs, fcn, is_haplotype_specific) {
  # convert sampled data to original data format
  # encode segment ID to keep the order (optional)
  cns_bs$sid <- 1:nrow(cns_bs)
  pos_names <- setdiff(names(cns_bs), c("sid"))
  
  if(is_haplotype_specific){
    cns_bs %>%
      gather(all_of(pos_names), key = "pos", value = "cn") %>%
      separate(pos, sep = "_", into = c("chr", "seg")) %>% rowwise() %>% 
      mutate(chr = as.integer(chr)) %>%
      arrange(sid, chr) -> dcn_bs0    
      # dcn_bs0 %>% arrange(sid, chr, seg) %>% View()
      # replace seg name so that the two haplotypes have the same name 
      # to ensure the same order, move "-2" to the end of the string
      dcn_bs0 %>% mutate(seg = ifelse(grepl("-2", seg), paste(str_remove(seg, "-2"), "-2", sep = ""), seg)) %>% arrange(sid, chr, seg) -> dcn_bs_rpl
      dcn_bs_rpl %>% mutate(seg = ifelse(grepl("-2", seg), str_replace(seg, "-2", ""), seg)) -> dcn_bs_same
      # check both haplotypes are there
      # dcn_bs_same %>% group_by(sid, chr, seg) %>% tally() %>% filter(n!=2)
      dcn_bs_same %>% group_by(sid, chr, seg) %>% mutate(seg = as.integer(seg), cn = paste0(cn, collapse = "-")) %>% separate(cn, into = c("cnA", "cnB"), sep = "-") -> dcn_bs
  }else{
    cns_bs %>%
      gather(all_of(pos_names), key = "pos", value = "cn") %>%
      separate(pos, sep = "_", into = c("chr", "seg")) %>%
      mutate(chr = as.integer(chr), seg = as.integer(seg)) %>%
      arrange(sid, chr) -> dcn_bs    
  }


  # check segment ID orders are the same across samples
  # dcn_bs %>% filter(sid == 1) %>% select(seg) -> seg1
  # dcn_bs %>% filter(sid == 2) %>% select(seg) -> seg2
  # all.equal(seg1, seg2)
  # dcn_bs %>% filter(sid == 1) %>% select(chr) -> chr1
  # dcn_bs %>% filter(sid == 2) %>% select(chr) -> chr2
  # all.equal(chr1, chr2)

  # write the data into gz file for tree building
  gz1 <- gzfile(fcn, "w")
  write.table(dcn_bs, gz1, quote = F, row.names = F, col.names = F, sep = "\t")
  close(gz1)
}


plot_tree <- function(tree) {
  p <- ggtree(tree)
  p <- p + geom_tiplab(align = TRUE) + theme_tree2()
  edge <- data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  colnames(edge) <- c("parent", "node", "edge_num", "edge_len")
  p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)
  return(p)
}


build_pratchet <- function(phydata) {
  # Build Maximum Parsimony tree
  mptree <- pratchet(phydata)
  mptree <- root(mptree, as.character(ncol(bp_df)), resolve.root = TRUE)
  mptree <- acctran(mptree, phydata)
  plot(mptree)
  # write.nexus(mptree, file=paste0("tree_mp", ".nhs"))
  write.tree(mptree, file = paste0("tree_mp", ".nwk"))
}


pick_top_ntree <- function(uniq_strees, num_select, dir_nwk, format = "NHS") {
  # Sort all trees by score
  scores <- c()
  for (i in 1:length(uniq_strees)) {
    # print(i)
    stree <- uniq_strees[[i]]
    # plot(stree)
    score <- attributes(stree)$pscore
    scores <- c(scores, score)
  }
  sindex <- order(scores)

  # Only take top trees
  num <- ifelse(length(uniq_strees) > num_select, num_select, length(uniq_strees))
  sel_trees <- list()
  for (i in 1:num) {
    j <- sindex[i]
    stree <- uniq_strees[[j]]
    sel_trees[[i]] <- stree
    # score = attributes(stree)$pscore
    if (format == "NWK") {
      write.tree(stree, file = file.path(dir_nwk, paste0("tree_mp_", i, ".nwk")))
    }
  }

  # all.equal.phylo(strees[[1]],strees[[3]])
  # comparePhylo(strees[[1]],strees[[3]])
  if (format == "NHS") {
    write.nexus(sel_trees, file = file.path(dir_nwk, "trees_mp.nhs"), translate = F)
  }
}


get_nwk_trees <- function(phydata, dir_nwk, num_generate, num_select, output_format) {
  strees <- list()
  for (i in seq(1, num_generate)) {
    stree <- random.addition(phydata, method = "fitch")
    stree <- root(stree, as.character(length(phydata)), resolve.root = TRUE)
    stree <- acctran(stree, phydata)
    # p = plot_tree(stree)
    # print(p)
    strees[[i]] <- stree
  }

  uniq_strees <- unique(strees)
  cat("There are", length(uniq_strees), "unique trees\n")

  dir.create(file.path(dir_nwk), showWarnings = FALSE)
  # Remove old files
  oldfiles <- dir(path = dir_nwk, pattern = "nwk")
  if (length(oldfiles) > 0) {
    file.remove(file.path(dir_nwk, oldfiles))
  }

  pick_top_ntree(uniq_strees, num_select, dir_nwk, output_format)
}


# Get input CNs, add normal sample if required
# for haplotype specific data, append data in another block 
get_cn <- function(file_cn, incl_normal, is_haplotype_specific) {
  data_cn <- read.table(file_cn)
  # only pick the first 4 or 5 columns
  if (is_haplotype_specific) {
    if (ncol(data_cn) > 5) {
      data_cn <- data_cn[, 1:5]
    }
    names(data_cn) <- c("sid", "chr", "seg", "cnA", "cnB")
    data_cn %>% select(sid, chr, seg, cn = cnA) -> cnA
    data_cn %>%
      select(sid, chr, seg, cn = cnB) %>%
      mutate(seg = paste0(seg, "-2")) -> cnB
    data_cn <- rbind(cnA, cnB)
  } else {
    if (ncol(data_cn) > 4) {
      data_cn <- data_cn[, 1:4]
    }
    names(data_cn) <- c("sid", "chr", "seg", "cn")
  }

  # data_cn %>% group_by(sid) %>% count()

  # Check if the last sample is normal
  nid <- length(unique(data_cn$sid))
  data_cn %>% filter(sid == nid) -> sn
  # ncn = unique(sn$cn)

  # if(length(ncn) == 1 && ncn[1] == 2){   // all CNs may be normal in some patient samples
  if (incl_normal) {
    cat("The last sample is normal\n")
  } else {
    # add a normal sample, use filter() to get segments
    data_cn %>% filter(sid == 1) -> s1
    if (is_haplotype_specific) {
      s1$cn <- 1
    } else {
      s1$cn <- 2
    }
    s1$sid <- length(unique(data_cn$sid)) + 1
    data_cn <- rbind(data_cn, s1)
  }

  return(data_cn)
}



option_list <- list(
  make_option(c("-f", "--file_cn"),
    type = "character", default = "",
    help = "input copy number file (only absolute copy numbers are allowed) [default=%default]", metavar = "character"
  ),
  make_option(c("-b", "--bootstrap"),
    type = "integer", default = 0,
    help = "Whether or not to generate bootstrap samples (0: not; 1: yes)  [default=%default]", metavar = "number"
  ),
  make_option(c("-i", "--input_format"),
    type = "integer", default = 0,
    help = "The input format for tree building (0: copy number; 1: breakpoint)  [default=%default]", metavar = "number"
  ),
  make_option(c("-o", "--output_format"),
    type = "integer", default = 0,
    help = "The output format of selected trees (0: NHS; 1: NWK)  [default=%default]", metavar = "number"
  ),
  make_option(c("-d", "--dir_nwk"),
    type = "character", default = "",
    help = "The directory to store results (NEWICK trees) [default=%default]", metavar = "character"
  ),
  make_option(c("-c", "--file_bs"),
    type = "character", default = "",
    help = "The file to store copy number data obtained from bootstrapping", metavar = "character"
  ),
  make_option(c("-m", "--incl_normal"),
    action = "store_true", default = FALSE,
    help = "Whether or not normal sample is included in the input [default=%default]"
  ),
  make_option(c("-a", "--is_haplotype_specific"), action = "store_true", default = FALSE, help = "Whether or not the input copy numbers are allele-specific [default=%default]"),
  make_option(c("-n", "--num_generate"),
    type = "integer", default = 100,
    help = "The number of trees to generate [default=%default]", metavar = "number"
  ),
  make_option(c("-s", "--num_select"),
    type = "integer", default = 100,
    help = "The maximum number of trees to select [default=%default]", metavar = "number"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
file_cn <- opt$file_cn
input_format <- opt$input_format
incl_normal <- opt$incl_normal
is_haplotype_specific <- opt$is_haplotype_specific
num_generate <- opt$num_generate
num_select <- opt$num_select
dir_nwk <- opt$dir_nwk
output_format <- ifelse(opt$output_format, "NWK", "NHS")
file_bs <- opt$file_bs
  
# input_format = 0
# num_generate = 100
# num_select = 100
# output_format = "NHS"
# incl_normal = F
# is_haplotype_specific = T
# original copy number data in the format of cnets
data_cn <- get_cn(file_cn, incl_normal, is_haplotype_specific)

if (opt$bootstrap) {
  cat("Generating bootstrapping data\n")
  data_bs <- get_bootstrap_phydata(data_cn, is_haplotype_specific)
  phydata <- data_bs$phydata
  cns_bs <- data_bs$cns_bs
  get_nwk_trees(phydata, dir_nwk, num_generate, num_select, output_format)
  get_bootstrap_cn(cns_bs, file_bs, is_haplotype_specific)
} else {
  if (input_format == 0) {
    cat("Using copy numbers to build tree\n")
    data_cn %>%
      unite(chr, seg, col = "pos", sep = "_") %>%
      spread(pos, cn) %>%
      select(-sid) -> cns_wide
    phydata <- get_phydata(cns_wide)
  } else {
    cat("Using breakpoints to build tree\n")
    phydata <- get_phydata_bp(data_cn)
  }

  get_nwk_trees(phydata, dir_nwk, num_generate, num_select, output_format)
}
