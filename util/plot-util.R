#!/usr/bin/env Rscript

# This file contains common functions to plot trees and copy numbers

suppressMessages(library(copynumber))
# suppressMessages(library(reshape))
suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(treeio))
suppressMessages(library(ape))
suppressMessages(library(tools))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
# suppressMessages(library(aplot))

# This script is used to plot phylogenetic trees:
# 1) plotting a single tree file with tips as 1 to n or with provided labels
# 2) plotting a single tree file with bootstrapping support
# 3) plotting all tree files in a directory

############### basic settings ###############
# cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0","#FCAE91", "#B9574E", "#76000D", "#8B0000", "#000000")
# Max CN to show in heatmap
MAX_CN = 6
# For absolute CN 0 to 8 (obtained from https://colorbrewer2.org, 8 classes, sequential data)
# cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0", '#fdd49e', '#fdbb84','#fc8d59','#ef6548','#d7301f','#990000')
# For absolute CN 0 to 6 (obtained from https://colorbrewer2.org, 5 classes, sequential data)
cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0", '#fdcc8a','#fc8d59','#e34a33','#b30000')
# For relative CN -4 to 4
cn_colors2 = c('#08519c','#3182bd', '#6baed6', '#9ecae1', "#f0f0f0", '#fdcc8a','#fc8d59','#e34a33','#b30000')

# xlab_cna_size = "Expected size of copy number alterations per site (Mbp)"
# xlab_cna_num = "Expected number of copy number alterations per site"
xlab_cna_age = "Patient age (years)"

MIN_BLEN = 1e-3
TIP_OFFSET = 0.5
BLEN_DIGIT = 2


############### functions to process input ###############
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


get.outfile.name <- function(tree_file, branch_num = 0){
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


get.cn.matrix <- function(cn_file){
  data_cn = read.table(cn_file)
  if(ncol(data_cn) == 4){
    names(data_cn) <- c("sid", "chr", "seg", "cn")
  }else{
    if(ncol(data_cn) != 5){
      stop("There must be either 4 or 5 columns!")
    }
    names(data_cn) <- c("sid", "chr", "seg", "cnA", "cnB")
    data_cn = data_cn %>% mutate(cn = cnA + cnB) %>% select(-c("cnA", "cnB"))
  }  
  data_cn %>% unite(chr, seg, col = "pos", sep = "_") %>% spread(pos, cn) %>% select(-sid) -> cns_wide
  
  return(cns_wide)
}


# add missing segments
get_normal_segs <- function(d_withpos, chr_end_arm, pos_file){
  samples = data.frame(sample = unique(d_withpos$sample))

  # data must be sorted by chr and start so that end = lead(start) works properly
  d_noindex = d_withpos %>% arrange(sample, chromosome, start) %>% select(-index)

  # get the end of intermediate intervals first
  df_invar = d_noindex %>% dplyr::mutate(end = lead(start) - 1)
  df_invar$start = d_withpos$end + 1
  df_invar = df_invar %>% filter(end >= start)
  if(nrow(df_invar) > 0){
    df_invar$cn = 2
  }

  # Update interval at the end of each chr
  # Find the current end in the data for each chr
  d_withpos %>% group_by(chromosome) %>% dplyr::summarise(start = max(end)) -> pos_max
  seg_max = merge(pos_max, chr_end_arm, by = c("chromosome")) %>% dplyr::mutate(start = start + 1) %>% select(chromosome, start, end = cend) %>% filter(start <= end)
  cnT_max = merge(seg_max, samples, all = T) %>% arrange(sample, chromosome, start)
  cnT_max$cn = 2

  # add intervals at the beginning of each chr
  d_withpos %>% group_by(chromosome) %>% dplyr::summarise(end = min(start)) -> pos_min
  seg_min = pos_min %>% dplyr::mutate(start = 1, end = end - 1) %>% filter(start <= end)
  cnT_min = merge(seg_min, samples, all = T) %>% arrange(sample, chromosome, start)
  cnT_min$cn = 2

  #�Remove interval at the end which has wrong boundaries
  max_ends = pos_max$start + 1
  min_ends = pos_min$end - 1
  row2del = df_invar %>% dplyr::filter(start %in% max_ends | end %in% min_ends)
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

  # d_noindex %>% group_by(chromosome, sample) %>% summarise(s = min(start), e = max(end)) %>% select(-sample) %>% unique() -> tmp
  # ftmp = str_replace(pos_file, "-cn.txt.gz", "-d_noindex.txt")
  # write_tsv(d_noindex, ftmp)
  #
  # if(nrow(cnT_min) > 0){
  # cnT_min %>% group_by(chromosome, sample) %>% summarise(s = min(start), e = max(end)) %>% select(-sample) %>% unique() -> tmp
  # ftmp = str_replace(pos_file, "-cn.txt.gz", "-cnT_min.txt")
  # write_tsv(cnT_min, ftmp)
  # }

  # df_invar_mid %>% group_by(chromosome, sample) %>% summarise(s = min(start), e = max(end)) %>% select(-sample) %>% unique() -> tmp
  # ftmp = str_replace(pos_file, "-cn.txt.gz", "-df_invar_mid.txt")
  # write_tsv(df_invar_mid, ftmp)

  # if(nrow(cnT_max) > 0){
  # cnT_max %>% group_by(chromosome, sample) %>% summarise(s = min(start), e = max(end)) %>% select(-sample) %>% unique() -> tmp
  # ftmp = str_replace(pos_file, "-cn.txt.gz", "-cnT_max.txt")
  # write_tsv(cnT_max, ftmp)
  # }
  #
  # if(nrow(seg_chm) > 0){
  # seg_chm %>% group_by(chromosome, sample) %>% summarise(s = min(start), e = max(end)) %>% select(-sample) %>% unique() -> tmp
  # ftmp = str_replace(pos_file, "-cn.txt.gz", "-seg_chm.txt")
  # write_tsv(seg_chm, ftmp)
  # }
  #
  # Add Seg for all regions
  cn_all = rbind(d_noindex, cnT_min, df_invar_mid, cnT_max, seg_chm) %>% arrange(sample, chromosome, start) %>% group_by(sample) %>% dplyr::mutate(index = row_number()) %>% select(sample, chromosome, index, start, end, cn)

  return(cn_all)
}


get.site.coord <- function(pos_file, cyto_file, bin_file = "", seed = NA){
  # original input, used to get sample names and original positions
  dpos <- read.table(pos_file)
  if(ncol(dpos)== 4 || ncol(dpos)== 5){  # the index has to start from 1 for each chromosome
    if(ncol(dpos)== 4){
      names(dpos) <- c("sample", "chromosome", "index", "cn")
    }else{
      names(dpos) <- c("sample", "chromosome", "index", "cnA", "cnB")
      dpos = dpos %>% mutate(cn = cnA + cnB) %>% select(-c("cnA", "cnB"))
    }
    
    nsite = dpos %>% group_by(sample) %>% tally() %>% select(n) %>% unique()
    if(bin_file != ""){   # read bin file to get positions
      if(file_ext(bin_file) ==  "Rdata"){
        bins = load(bin_file)
        if(nsite == 4401){
          chr_sites_full = bins_4401 
        }else{
          data = dpos %>% spread(sample, cn)
          bin_count = data %>% group_by(chromosome) %>% count() 
          nested_bins <- bins_4401 %>%
            group_by(chromosome) %>%   # prep for work by Species
            nest() %>%              # --> one row per Species
            ungroup() %>%
            mutate(n = bin_count$n) # add sample sizes
          
          sampled_bins <- nested_bins %>%
            mutate(samp = map2(data, n, sample_n))
          
          chr_sites_full = sampled_bins %>%
            select(-data) %>%
            unnest(samp) %>% arrange(chromosome, start) %>% select(-c(n))
        }
      }else{
        chr_sites_full = readRDS(bin_file)
        names(chr_sites_full) = c("chromosome", "start", "end")
      }
    }else{
      # randomly split the genome into n sites
      chr_end_arm = get_hg_ends(cyto_file)
      # get number of sites for each chromosome
      dpos %>% filter(sample == 1) %>% group_by(chromosome) %>% tally() -> nsite_per_chr
      chr_sites = merge(chr_end_arm, nsite_per_chr, by = c("chromosome"))
      chr_sites = chr_sites %>% select(chromosome, cend, n)

      if(!is.na(seed)) set.seed(seed)
      chr_sites_all = chr_sites %>% rowwise() %>% mutate(sites = paste(sample.int(cend, n), collapse  = ","))

      max_nsite = max(chr_sites$n)
      snames = paste0("s", 1:max_nsite)
      chr_sites_sep = chr_sites_all %>% separate(sites, sep = ",", into = snames)
      chr_site_long = chr_sites_sep %>% gather(all_of(snames), key = "boundary", value = "start") %>% drop_na() %>% select(-boundary) %>% mutate(start = as.integer(start)) %>% arrange(chromosome, start)
      chr_site_long$end = lead(chr_site_long$start) - 1
      # adjust the boundary of first and last boundary
      site_ends = chr_site_long %>% group_by(chromosome) %>% summarise(max_end = max(start), min_end = min(start))
      chr_sites_full = merge(chr_site_long, site_ends, by = c("chromosome"))
      chr_sites_full = chr_sites_full %>% rowwise() %>%  mutate(start = ifelse(start == min_end, 1, start), end = ifelse(start == max_end, cend, end))
    }

    # extract position and merge with input file
    chr_sites_pos = chr_sites_full  %>% group_by(chromosome) %>% mutate(index = 1:n()) %>% select(chromosome, index, start, end)
    dpos = merge(dpos, chr_sites_pos, by = c("chromosome", "index"))
    dpos = dpos %>% select(sample, chromosome, index, cn, start, end) %>% arrange(sample, chromosome, index)

  }else{
    names(dpos) <- c("sample", "chromosome", "index", "cn", "start", "end")
  }

  return(dpos)
}


# convert copy number state to haplotype-specific copy number
state2hcn <- function(state, cn_max = 4){
  sums = c()
  nmax = cn_max + 1
  for(i in 1:nmax){
    s = i * (i + 1) / 2
    sums = c(sums, s)
  }
  
  if(state < sums[1]) return(c(0, 0));
  
  i = 2;
  repeat{
    if(state >= sums[i - 1] & state < sums[i]){
      cnA = state - sums[i - 1];
      cnB = i - 1 - cnA;
    }
    i = i + 1;
    if(i > nmax) break;
  }
  
  
  return(c(cnA, cnB))
}


# convert copy number state to a total copy number
state2tcn <- function(state, cn_max = 4){
  sums = c()
  nmax = cn_max + 1
  for(i in 1:nmax){
    s = i * (i + 1) / 2
    sums = c(sums, s)
  }
  
  if(state < sums[1]) return(0);
  
  i = 2;
  repeat{
    if(state >= sums[i - 1] & state < sums[i]){
      cnA = state - sums[i - 1];
      cnB = i - 1 - cnA;
    }
    i = i + 1;
    if(i > nmax) break;
  }
  
  tcn = cnA + cnB
  
  return(tcn)
}

# for(i in 0:9){
#   print(i)
#   # print(state2tcn((i)))
#   print(state2hcn((i)))
# }


# has_normal: whether or not the normal sample is in the input
get.all.state <- function(dans, cn_file, seg_file, dpos, pos_file, has_normal = T, is_haplotype_specific = F, cn_max = 4){
  dpos %>% select(chromosome, index, start, end) %>% unique() -> site_pos

  if(seg_file != ""){  # reconstructed states
    samples = unique(dpos$sample) %>% unlist()

    dseg <- read.table(seg_file)
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
    if(has_normal) samples = samples[1 : (length(samples) - 1)]
    names(dseg_cn) = c("chromosome", "index", paste0("s", all_of(samples)))

    dseg_cn_long = dseg_cn %>% gather(paste0("s", all_of(samples)), key = "sample", value = "cn") %>% mutate(sample = str_replace(sample, "s","")) %>% group_by(sample, chromosome)  %>%  mutate(index = 1:n()) %>% select(sample, chromosome, index, cn)
    dseg_cn_long$sample = as.numeric(dseg_cn_long$sample)
    
    if(is_haplotype_specific){
      dseg_cn_long = dseg_cn_long %>% rowwise() %>% mutate(cn = state2tcn(cn, cn_max))
    }
    # dseg_cn_long %>% group_by(sample) %>% tally()

    # add normal samples as seg_file does not contain normal samples
    dseg_cn_long %>% filter(sample == 1) -> s1  # use filter() to get segments
    s1$cn = 2
    s1$sample = length(unique(dseg_cn_long$sample)) + 1
    dtip = rbind(dseg_cn_long, s1)

    dstate = rbind(dans, dtip)
    # get reconstructed states with position
    d_all = merge(dstate, segs_pos_all, by = c("chromosome", "index"), all.x = T) %>% select(all_of(names(dpos))) %>% arrange(sample, chromosome, index)
    # d_all %>% filter(sample == 2 & chromosome == 1 & cn != 2) %>% print()

  }else{ # simulated states (normal sample is included by cnets)
    if(pos_file == ""){
      stop("either segment file or input file must be provided!")
    }
    dtip = read.table(pos_file)
    if(ncol(dtip) > 4) dtip = dtip[,1:4]
    names(dtip) <- c("sample", "chromosome", "index", "cn")

    if(!has_normal){
      dtip %>% filter(sample == 1) -> s1  # use filter() to get segments
      s1$cn = 2
      s1$sample = length(unique(dtip$sample)) + 1
      dtip = rbind(dtip, s1)
    }

    ntip = length(unique(dtip$sample))
    root = ntip + 1
    dans = dans %>% filter(sample != root)

    dstate = rbind(dans, dtip)
    # get reconstructed states with position
    d_all = merge(dstate, site_pos, by = c("chromosome", "index"), all.x = T) %>% select(all_of(names(dpos))) %>% arrange(sample, chromosome, index)
    # d_all %>% filter(sample == 2 & chromosome == 1 & cn != 2) %>% print()
  }

  return(d_all)
}


# combine CNP with position information to get a file with the format: sample, chrom, start, end, cn
# use chr and site index to bind the two datasets
# ref_file contains the reference position of each site
get.cn.data.by.pos <- function(cn_file, pos_file, seg_file, cyto_file, labels, ordered_nodes, has_normal = F, bin_file = "", seed = NA, is_haplotype_specific = F, cn_max = 4, excluded_tip = ""){
  dans <- read.table(cn_file)
  if(ncol(dans) == 4){
    names(dans) <- c("sample", "chromosome", "index", "cn")
  }else{
    if(ncol(dans) != 5){
      stop("There must be either 4 or 5 columns!")
    }
    names(dans) <- c("sample", "chromosome", "index", "cnA", "cnB")
    dans = dans %>% mutate(cn = cnA + cnB) %>% select(-c("cnA", "cnB"))
  }

  npos_ans = dans %>% group_by(sample) %>% tally() %>% select(n) %>% unique()

  # extract positions for all sites
  dpos = get.site.coord(pos_file, cyto_file, bin_file, seed)
  # print(dpos)
  # print(summary(dpos))
  # number of sites in original input
  npos = dpos %>% group_by(sample) %>% tally() %>% select(n) %>%  unique()

  # get reconstructed states with position
  # read segment file for original input as invariable bins are excluded and consecutive bins may be merged
  d_all = get.all.state(dans, cn_file, seg_file, dpos, pos_file, has_normal, is_haplotype_specific, cn_max)
  # d_all %>% filter(sample %in% unique(dpos$sample)) %>% summary() %>% print()

  # d_all %>% group_by(chromosome, sample) %>% tally() %>% select(-sample) %>% unique() %>% print()
  # d_all %>% group_by(chromosome, sample) %>% summarise(s = min(start), e = max(end)) %>% select(-sample) %>% unique() -> tmp
  # ftmp = str_replace(pos_file, "-cn.txt.gz", "-chr-icmpl.txt")
  # write_tsv(tmp, ftmp)

  # add missing normal segments
  chr_end_arm = get_hg_ends(cyto_file)
  # ftmp = str_replace(pos_file, "-cn.txt.gz", "-hg19.txt")
  # write_tsv(chr_end_arm, ftmp)
  # print(d_all)
  if(npos_ans < npos){
    d_all = get_normal_segs(d_all, chr_end_arm, pos_file)
  }
  # print("after getting normal segments")
  # print(d_all)
  # replace node IDs with sample names
  ans_nodes = unique(dans$sample) %>% unlist()
  ans_labels = data.frame(sample = ans_nodes, name = ans_nodes)
  
  if(length(labels) > 0){
    tip_labels = data.frame(sample = 1:length(labels), name = labels)
  }else{
    if(has_normal){
      tip_labels = data.frame(sample = 1:length(unique(dpos$sample)), name = unique(dpos$sample)) 
    }else{
      ns = length(unique(dpos$sample)) + 1
      tip_labels = data.frame(sample = 1:ns, name = c(unique(dpos$sample), ns)) 
    }
  }
  if(excluded_tip != ""){
    tip_labels = tip_labels %>% filter(!name %in% excluded_tip)
  }
  # print(tip_labels)
  node_labels = rbind(ans_labels, tip_labels)
  
  d_all = merge(d_all, node_labels, by = c("sample"), all.y = T) %>% select(-sample) %>% select(sample = name, chromosome, index, start, end, cn)
  d_all$sample = factor(d_all$sample, levels = ordered_nodes)

  # d_all %>% group_by(chromosome , sample) %>% tally() %>% select(-sample) %>% unique() %>% print()
  # d_all %>% group_by(chromosome, sample) %>% summarise(s = min(start), e = max(end)) %>% select(-sample) %>% unique() -> tmp
  # ftmp = str_replace(pos_file, "-cn.txt.gz", "-chr.txt")
  # write_tsv(tmp, ftmp)

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
get.bootstrap.tree <- function(mytree, labels, bstrap_dir, pattern, excluded_tip = ""){
  btrees = list()
  cat("Patterns to match bootstrapping trees: ", pattern, "\n")
  files = list.files(path = bstrap_dir, pattern = glob2rx(pattern), recursive = F)
  for(i in 1:length(files)){
    fname = file.path(bstrap_dir, files[i])
    if(endsWith(fname, ".txt")){
      dt = read.table(fname,header = T)
      btree = make.tree(dt, labels)    
    }else{
      if(!endsWith(fname, ".nex")){
        stop("The tree files should be in nexus format, with file name ending with .nex!")
      }
      btree = read.nexus(fname)
      if(length(labels) > 0){
        lbl_orders = 1:length(labels)
        btree$tip.label = labels[match(btree$tip.label, lbl_orders)]        
      }
    }
    
    if(excluded_tip != ""){
      btree = drop.tip(btree, excluded_tip)
    } 
    
    btrees[[i]] = btree
  }

  # prop.clades calls internally prop.part with the option check.labels = TRUE
  # prop.part counts the number of bipartitions found in a series of trees given as .... If a single tree is passed, the returned object is a list of vectors with the tips descending from each node (i.e., clade compositions indexed by node number).
  clad <- prop.clades(mytree, btrees, rooted = TRUE)
  clad[is.na(clad)] = 0
  bsval = as.integer(clad * 100 / length(files))
  if(length(mytree$node.label) > 0){
    mytree$node.label = paste0(mytree$node.label, "[&bootstrap=", bsval, "]")
  }else{
    mytree$node.label = bsval
  }
 
  return(mytree)
}


# get the tree with confidence intervals at internal nodes
# need to ensure both the original tree and bootstrap trees are read in the same way
# when reading txt files, the internal node IDs are not recorded
# when reading nex files, the internal node IDs are different from those in the file
# tree_file_nex: bootstrap tree or original tree in nexus format
# tip.label: a vector of mode character giving the labels of the tips; the order of these labels corresponds to the integers 1 to n in edge.
# node.label (optional) a vector of mode character giving the labels of the nodes (ordered in the same way than tip.label).
get.ci.tree <- function(tree_file_nex, bstrap_dir, labels, has_bstrap = F, nex_pattern = "*.nex", ci_prefix = "time_0.95_CI", excluded_tip = ""){
  bstrees = list.files(bstrap_dir, pattern = nex_pattern)
  ntimes_all = data.frame()   # not always time, depending on the meaning of branch lengths
  nbs = length(bstrees)

  for(i in 1:nbs){
    # i = 2
    ftree = file.path(bstrap_dir, bstrees[i])
    btree = read.nexus(ftree)
    if(excluded_tip != ""){
      btree = drop.tip(btree, excluded_tip)
    }     
    # print(btree)
    # get node time
    ndepths = node.depth.edgelength(btree)
    max_depth = max(ndepths)
    ntimes = as.data.frame(max_depth - ndepths)
    colnames(ntimes) = ""
    ntimes_all = rbind(ntimes_all, t(ntimes))
  }
  # print(ntimes_all)
  # map Id to node label
  # all the bootstrapping trees have the same topology here
  id_labels = c(btree$tip.label, btree$node.label)
  colnames(ntimes_all) = id_labels

  # if(excluded_tip != ""){
  #   stree = read.nexus(tree_file_nex)
  #   stree = drop.tip(stree, excluded_tip)
  #   fout = str_replace(tree_file_nex, ".nex", "_dtip.nex")
  #   write.nexus(stree, file=fout)
  #   stree = read_file(fout)
  #   mega_tree = read.mega(fout)
  # }else{
    # bootstrap tree already has relevant tip removed in the previous step
    # read the tree string
    stree = read_file(tree_file_nex)  # both original tree and bootstrap tree have the same node labels
    # keep bootstrap value with read.mega
    mega_tree = read.mega(tree_file_nex)    
  # }

  # compute CI interval and append CI to the tree for visualization
  start_inode = btree$Nnode + 2
  root = start_inode

  for(i in start_inode: ncol(ntimes_all)){  # id increasing from root
    # i = 8
    lbl = id_labels[i]   # same labels as the original tree
    # print(i)
    # print(lbl)
    times = ntimes_all[, lbl]

    # smu = mean(times)
    # error = qnorm(0.975) * sd(times) / sqrt(length(times))
    # left = smu - error
    # right = smu + error
    sorted_times = sort(times)
    # print(sorted_times)
    ci95 = quantile(sorted_times, probs = c(0.025, 0.975))
    left = ci95[1]
    right = ci95[2]

    # when the tree has bootstrap support value, the node label is the support value, assume: all internal branches have different lengths, not work at extreme cases such as when showing branch length as number of mutations and many branches have none
    # [&time_95%_CI={4.50612194676725E-002,1.10259531669597E-001}]
    if(i == start_inode){
      marker = ";"
      # blen = ""
    }else{
      marker = ":"
      # use branch length to avoid mismatching
      # find the edge ID with the node
      # eid = which(mega_tree@phylo$edge[, 2]==i)
      # blen = mega_tree@phylo$edge.length[eid]
    }

    if(has_bstrap){
      # find the index of current node
      # internal node are relabeled,
      lid = which(mega_tree@phylo$node.label == lbl)
      bid = which(mega_tree@data$node == lid + mega_tree@phylo$Nnode + 1)  # mega_tree@data$node are node IDs
      bsval = mega_tree@data$bootstrap[bid]
      orig_str = paste0(lbl, "\\[&bootstrap=", bsval, "]", marker)
      ci_str = paste0(lbl, "\\[&", ci_prefix, "={", left, ",", right, "},bootstrap=", bsval, "]", marker)
    }else{
      orig_str = paste0(lbl, marker)
      ci_str = paste0(lbl, "[&", ci_prefix, "={", left, ",", right, "}]", marker)
    }
    # print(orig_str)
    # print(ci_str)
    stree = str_replace(stree, orig_str, ci_str)
  }

  # write annotated tree to a file
  bname = basename(tree_file_nex)
  dname = dirname(tree_file_nex)
  prefix = strsplit(bname, "\\.")[[1]][1]
  ftree_ci = file.path(dname, paste0(prefix, ".ci.nex"))
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

  # For BE data
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


###########################  functions to plot trees ###########################

plot.tree <- function(tree, title = "", da = data.frame(), add_blen = T) {
  p <- ggtree(tree)  #+ geom_rootedge()
  p <- p + geom_tiplab()
  p <- p + geom_text2(aes(subset = !isTip, label = label), hjust = -0.3)
  p <- p + ggtitle(title)

  if(add_blen){
    edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
    colnames(edge) = c("parent", "node", "edge_num", "edge_len")
    edge$edge_len = round(edge$edge_len, BLEN_DIGIT)
    p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)
  }

  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }

  return(p)
}


# Plot tree with xlim specified to show full tip labels (representing time backwards)
plot.tree.xlim <- function(tree, title = "", rextra = 20, da = data.frame(), add_blen = T) {
  p <- ggtree(tree, size = 0.5, linetype = 1)
  # Add margin to show full name of labels if (is.na(tree.max))
  tree.max = max(node.depth.edgelength(tree)) + rextra
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + geom_treescale() + xlim(NA, tree.max)
  p <- p + ggtitle(title)

  if(add_blen){
    edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
    colnames(edge) = c("parent", "node", "edge_num", "edge_len")
    edge$edge_len = round(edge$edge_len, BLEN_DIGIT)
    p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)
  }

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
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + theme_tree2() + xlim(xl, tree_max) + xlab(xlab_cna_age)
  p <- p + geom_text2(aes(subset=!isTip, label=label), hjust=-.3)

  if(nrow(da) > 0){
    p = get.plot.annot(tree, da, p)
  }

  return(p)
}


plot.tree.bootstrap <- function(tree, title = "", rextra = 20, da = data.frame(), scale_factor = 1){
  p <- ggtree(tree)
  #+ geom_rootedge()

  # support <- character(length(tree$node.label))
  # #The following three lines define your labeling scheme.
  # support[tree$node.label >= 95] <- "red"
  # support[tree$node.label < 95 & tree$node.label >= 70] <- "pink"
  # support[tree$node.label < 70] <- "blue"

  tree.max = max(node.depth.edgelength(tree)) + rextra
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + geom_treescale() + xlim(NA, tree.max)
  p <- p + geom_text2(aes(subset=!isTip, label=label), hjust=-.3, color="red")

  # edge = data.frame(tree$edge, edge_num = 1:nrow(tree$edge), edge_len = tree$edge.length)
  # colnames(edge)=c("parent", "node", "edge_num", "edge_len")
  # p <- p %<+% edge + geom_text(aes(x = branch, label = edge_len), nudge_y = 0.1)
  p + ggtitle(title)

  # if(scale_factor == 1){
  #   p <- p + xlab(xlab_cna_num)
  # }else{
  #   p <- p + xlab(xlab_cna_size)
  # }

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
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + theme_tree2() + xlim(xl, tree_max) + xlab(xlab_cna_age)
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
  # use patient age to shift x-axis
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

  p <- p + geom_range('time_0.95_CI', color = 'red', size = 3, alpha = 0.3)
  p <- p + xlim(xl, tree_max) + xlab(xlab_cna_age)

  if(nrow(da) > 0){
    p = get.plot.annot(tree_ci@phylo, da, p)
  }

  return(p)
}


# Plot CI for mutation number (size) estimation
plot.tree.ci.node.mut <- function(tree_ci, title = "", lextra = 0, rextra = 3, da = data.frame(), has_bstrap = F, has_inode_label = F, scale_factor = 1){
  tree_depth = max(node.depth.edgelength(tree_ci@phylo))
  cat("tree depth ", tree_depth, "\n")
  tree_max = tree_depth + rextra
  xl = 0 - lextra

  p <- ggtree(tree_ci)
  p <- p + geom_tiplab(align = TRUE, offset = TIP_OFFSET) + geom_treescale() + ggtitle(title)

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

  if(scale_factor == 1){
    p <- p + geom_range('nmut_0.95_CI', color='red', size = 3, alpha = 0.3)
    p <- p + xlim(xl, tree_max) 
    # + xlab(xlab_cna_num)
  }else{
    p <- p + geom_range('mutsize_0.95_CI', color='red', size=3, alpha = 0.3)
    p <- p  + xlim(xl, tree_max) 
    # + xlab(xlab_cna_size)
  }


  if(nrow(da) > 0){
    p = get.plot.annot(tree_ci@phylo, da, p)
  }

  return(p)
}


# Plot CI for mutation number (size) estimation with a simplified layout: no tip alignment
plot.tree.ci.node.mut.smpl <- function(tree_ci, title = "", lextra = 0, rextra = 3, da = data.frame(), has_bstrap = F, has_inode_label = F, scale_factor = 1){
  tree_depth = max(node.depth.edgelength(tree_ci@phylo))
  tree_max = tree_depth + rextra
  xl = 0 - lextra

  p <- ggtree(tree_ci)
  p <- p + geom_tiplab() + ggtitle(title) + geom_treescale()

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

  if(scale_factor == 1){
    p <- p + geom_range('nmut_0.95_CI', color='red', size = 3, alpha = 0.3)
    p <- p + xlim(xl, tree_max)
  }else{
    p <- p + geom_range('mutsize_0.95_CI', color='red', size=3, alpha = 0.3)
    p <- p + xlim(xl, tree_max)
  }


  if(nrow(da) > 0){
    p = get.plot.annot(tree_ci@phylo, da, p)
  }

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


###########################  functions to plot copy numbers ###########################
# Plot CNPs of multiple samples as a heatmap for easy comparison
# Format of d_seg: sample, chrom, start, end, cn, ploidy
plot.cn.heatmap <- function(d_seg, main, type="absolute", theme = theme1, cn_colors = cn_colors1, allele_specific = F){
  d_seg$cn = round(d_seg$cn)
  d_seg$chrom = factor(d_seg$chrom, levels=paste("",c(c(1:22, "X")),sep=""))
  d_seg$pos = (d_seg$end + d_seg$start) / 2
  d_seg$width = (d_seg$pos - d_seg$start) * 2

  if(type=="absolute"){
    print("Plot absolute copy number")
    d_seg %>% dplyr::mutate(cn=if_else(cn > MAX_CN, MAX_CN, cn)) -> d_seg
    cn_vals = c("0", "1", "2", "3", "4", "5", "6")
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



get.cn.data.by.bin <- function(d, bins_4401){
  data = d %>% spread(sample, cn)
  # md <- melt(d, id=c("sample","chromosome","index"))
  # data <- cast(md, chromosome + index ~ sample)

  data %>% group_by(chromosome) %>% count() -> bin_count
  bins_4401 %>% group_by(chromosome) %>% count() -> ref_bin_count

  if(sum(bin_count$n) == sum(ref_bin_count$n)){
    # When the input data have 4401 bins
    data <- cbind(bins_4401, data)
  }
  else{ # When the data have less bins, randomly pick these bins among all for simplicity, not covering whole genome
    # Change reference bins to the desired number
    # adapted from https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html
    nested_bins <- bins_4401 %>%
      group_by(chromosome) %>%   # prep for work by Species
      nest() %>%              # --> one row per Species
      ungroup() %>%
      mutate(n = bin_count$n) # add sample sizes

    sampled_bins <- nested_bins %>%
      mutate(samp = map2(data, n, sample_n))

    sampled_bins %<>%
      select(-data) %>%
      unnest(samp) %>% arrange(chromosome, start) %>% select(-c(n))

    data <- cbind(sampled_bins, data)
  }

  data <- data[,-c(3,4,5)]

  return(data)
}


# Plot CNP as segment lines
plot.cn <- function(in_file, out_file, bins_4401){
  d <- read.table(in_file, header=FALSE)
  names(d) <- c("sample", "chromosome", "index", "cn")
  nsample <- length(unique(d$sample))

  data <- get.cn.data.by.bin(d, bins_4401)

  #par(ask=F)
  pdf(out_file, height=10, width=20)
  plotGenome(data=data, ylab="copy number", layout=c(2,3), sample=c(1:nsample), equalRange=FALSE, q=0, col="red" )
  dev.off()
}
