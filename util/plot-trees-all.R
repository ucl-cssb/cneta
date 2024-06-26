#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

# from https://stackoverflow.com/questions/47044068/get-the-path-of-current-script/47045368
getCurrentFileLocation <-  function(){
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

source(file.path(getCurrentFileLocation(), "plot-util.R"))


option_list = list(
  make_option(c("-f", "--tree_file"), type="character", default="",
              help="tree file name in TSV format [default=%default]", metavar="character"),
  make_option(c("", "--tree_file_nex"), type="character", default="",
              help="tree file name in NEXUS format [default=%default]", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default="",
              help="The name of output file [default=%default]", metavar="character"),
  make_option(c("-a", "--annot_file"), type="character", default="",
              help="The file containing the labels of tip nodes with at least two columns with header name 'sample'(or 'Sample') and 'id', where ids are ordered from 1 to the number of samples [default=%default]", metavar="character"),
  make_option(c("", "--time_file"), type="character", default="",
              help="The file containing the sampling time information [default=%default]", metavar="character"),
  make_option(c("", "--cn_file"), type="character", default="",
              help="The file which contains the copy numbers of internal nodes [default=%default]", metavar="character"),
  make_option(c("", "--pos_file"), type="character", default="",
              help="The file which contains the positions of each site along with the copy numbers of tip nodes (the index needs to start from 1 for each chromosome) [default=%default]", metavar="character"),
  make_option(c("", "--seg_file"), type="character", default="",
              help="The file which contains the segment start and end ID in the input copy number file [default=%default]", metavar="character"),
  make_option(c("", "--cyto_file"), type="character", default="",
              help="The file which contains the chromosome boundaries in human reference genome (e.g. hg19) [default=%default]", metavar="character"),
  make_option(c("", "--bin_file"), type="character", default="",
              help="The RDS file which contains the positions of bins used for calling copy numbers in human reference genome (e.g. hg19) with 3 columns: chr, start, end [default=%default]", metavar="character"),
  make_option(c("-d", "--tree_dir"), type="character", default="",
              help="The directory containing all the tree files to plot [default=%default]", metavar="character"),
  make_option(c("-s", "--bstrap_dir"), type="character", default="",
              help="The directory containing all the bootstrapping tree files for tree inference [default=%default]", metavar="character"),
  make_option(c("", "--bstrap_dir2"), type="character", default="",
              help="The directory containing all the bootstrapping tree files when tree topology is fixed [default=%default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="",
              help="The naming pattern of tree files [default=%default]", metavar="character"),
  make_option(c("", "--nex_pattern"), type="character", default=".nex",
              help="The naming pattern of bootstrap tree files when topology is fixed, used to get node time confidence intervals [default=%default]", metavar="character"),
	make_option(c("-b", "--branch_num"), type="integer", default = 0,
              help="The type of values on branches (0: time in year, 1: mutation number) [default=%default]", metavar="integer"),
  make_option(c("", "--scale_factor"), type="numeric", default = 1,
              help="The scale factor to convert branch length unit from #mutations to size of copy number changes [default=%default]", metavar="numeric"),
  make_option(c("", "--lextra"), type="numeric", default = 3,
              help="The left margin for the plot [default=%default]", metavar="numeric"),
  make_option(c("", "--rextra"), type="numeric", default = 3,
              help="The left margin for the plot [default=%default]", metavar="numeric"),
  make_option(c("", "--width"), type="numeric", default = 8,
              help="The width for the plot [default=%default]", metavar="numeric"),              
  make_option(c("", "--height"), type="numeric", default = 5,
              help="The height for the plot [default=%default]", metavar="numeric"),
  make_option(c("", "--font_size"), type="numeric", default = 6,
              help="The font size for labels on the heatmap plot [default=%default]", metavar="numeric"),              
  make_option(c("-m", "--mut_rate"), type="numeric", default = 0,
              help="The rate of somatic chromosomal aberrations [default=%default]", metavar="numeric"),
  make_option(c("-u", "--dup_rate"), type="numeric", default = 0,
              help="The site duplication rate of somatic chromosomal aberrations [default=%default]", metavar="numeric"),
  make_option(c("-e", "--del_rate"), type="numeric", default = 0,
              help="The site deletion rate of somatic chromosomal aberrations [default=%default]", metavar="numeric"),
  make_option(c("", "--is_haplotype_specific"), default = FALSE, action = "store_true",
              help="The input copy number is haplotype-specific [default=%default]", metavar="logical"),
  make_option(c("", "--cn_max"), type="integer", default = 4,
              help="The maximum copy number specified in building the tree, used for plotting haplotype-specific input [default=%default]", metavar="integer"),  
  make_option(c("-w", "--with_title"), default = FALSE, action = "store_true",
              help="Showing title [default=%default]", metavar="logical"),
  make_option(c("", "--with_cn"), default = FALSE, action = "store_true",
              help="Plotting copy numbers [default=%default]", metavar="logical"),
  make_option(c("", "--has_normal"), default = FALSE, action = "store_true",
              help="Whether the input has copy numbers for normal sample [default=%default]", metavar="logical"),
  make_option(c("", "--excluded_tip"), type="character", default="",
              help="Tips to exclude when plotting the tree, string separated by ',' if there are more than one tip [default=%default]", metavar="character"),
  make_option(c("", "--title"), type="character", default="",
              help="The title of the plot [default=%default]", metavar="character"),
  make_option(c("-t", "--plot_type"), type="character", default="single",
              help="The type of plot, including: all (plotting all tree files in a directory), single (plotting a single tree file), bootstrap (plotting a single tree file with bootstrapping support) [default=%default]", metavar="character"),
  make_option(c("-l", "--tree_style"), type="character", default="simple",
              help="The style of tree plot, including: simple (a simple tree with tip labels and branch lengths), xlim (adding xlim to the tree), age (x-axis as real age of the patient), and ci (plotting a single tree file with confidence interval of node ages) [default=%default]", metavar="character"),
  make_option(c("", "--ggtree_style"), type="integer", default = 0,
              help="The style of showing ggtree (0: default, 1: simplified (used for plot with bootstrap support values)) [default=%default]", metavar="integer"),
  make_option(c("", "--seed"), type="numeric", default = NA,
              help="The seed used for sampling sites in the genomes [default=%default]", metavar="numeric")

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
bin_file = opt$bin_file
plot_type = opt$plot_type
tree_dir = opt$tree_dir
pattern = opt$pattern
nex_pattern = opt$nex_pattern
bstrap_dir = opt$bstrap_dir
bstrap_dir2 = opt$bstrap_dir2
tree_style = opt$tree_style
ggtree_style = opt$ggtree_style
annot_file = opt$annot_file
branch_num = opt$branch_num
scale_factor = opt$scale_factor
mut_rate = opt$mut_rate
dup_rate = opt$dup_rate
del_rate = opt$del_rate
is_haplotype_specific = opt$is_haplotype_specific
cn_max = opt$cn_max
with_title = opt$with_title
with_cn = opt$with_cn
has_normal = opt$has_normal
title = opt$title
seed = opt$seed
height = opt$height
width = opt$width
font_size = opt$font_size
excluded_tip = str_split(opt$excluded_tip, ",")[[1]]

# customized font size 
theme1 = theme(legend.position = "none",
               #strip.text.x = element_blank(),
               #strip.text.y = element_blank(),
               strip.text.y.left = element_text(size=font_size, angle=0),
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

if(!is.na(seed)) set.seed(seed)

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
    out_file = get.outfile.name(fname, branch_num)
    mytree = get.tree(fname, branch_num, labels, scale_factor)
    if(excluded_tip != ""){
      mytree = drop.tip(mytree, excluded_tip)
    }
    p = print.single.tree(mytree, tree_style, time_file, title, lextra, rextra, da)
    ggsave(out_file, p, width = width, height = height)
  }

}else if(plot_type == "single"){
  cat(paste0("Plotting the tree in ", tree_file, "\n"))
  if(out_file == ""){
    out_file = get.outfile.name(tree_file, branch_num)
  }

  if(tree_style != "ci"){
    if(file_ext(tree_file) == "txt"){
      mytree = get.tree(tree_file, branch_num, labels, scale_factor)
    }else{
      mytree = read.nexus(tree_file)
    }
    
    if(excluded_tip != ""){
      mytree = drop.tip(mytree, excluded_tip)
    }    
    p = print.single.tree(mytree, tree_style, time_file, title, lextra, rextra, da)
  }else{
    # cat(paste0("Plotting confidence intervals for the tree in ", tree_file_nex, " with bootstrap trees in folder ", bstrap_dir2, "\n"))
    has_bstrap = F # assume no branch support value in the tree
    if(branch_num == 0){
      tree_ci = get.ci.tree(tree_file_nex, bstrap_dir2, labels, has_bstrap, nex_pattern, excluded_tip = excluded_tip)
      p = plot.tree.ci.node(tree_ci, time_file, title, lextra, rextra, da, T, T)
    }else{
      if(scale_factor == 1){
        ci_prefix = "nmut_0.95_CI"
      }else{
        ci_prefix = "mutsize_0.95_CI"
      }
      mytree = get.ci.tree(tree_file_nex, bstrap_dir2, labels, has_bstrap, nex_pattern, ci_prefix, excluded_tip = excluded_tip)
      p = plot.tree.ci.node.mut(mytree, time_file, title, lextra, rextra, da, T, T, scale_factor)
    }
  }
  
  if(with_cn){
      d = fortify(mytree)
      ordered_nodes = d$label[order(d$y, decreasing = T)]
      d_seg = get.cn.data.by.pos(cn_file, pos_file, seg_file, cyto_file, labels, ordered_nodes, has_normal, bin_file, seed, is_haplotype_specific, cn_max, excluded_tip)
      # get the node order of the tree and reorder heatmap
      d_seg = d_seg %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))
      phmap = plot.cn.heatmap(d_seg, "")

      pc = ggarrange(p, phmap, nrow = 1, widths = c(6.5, 9.5))
      # pc = phmap %>% insert_left(p, width = 0.6)  # not work for internal nodes
      #  width = 16
      ggsave(out_file, pc, width = width, height = height)
  }else{
      ggsave(out_file, p, width = width, height = height)
  }

}else if(plot_type == "bootstrap"){
  cat(paste0("Plotting bootstrap values for the tree in ", tree_file, "\n"))

  # when the tree is read from txt file, the internal node labels follow the same order as input CNs, but not explicit
  # mytree = get.tree(tree_file, branch_num, labels, scale_factor)
  mytree = read.nexus(tree_file_nex)
  if(length(labels) > 0){
    lbl_orders = 1:length(labels)
    mytree$tip.label = labels[match(mytree$tip.label, lbl_orders)]   
  }
  if(excluded_tip != ""){
    mytree = drop.tip(mytree, excluded_tip)
  }   
  # need to ensure the type of the branches in ML tree and bootstrap trees are the same (either time or number of mutations)
  mytree = get.bootstrap.tree(mytree, labels, bstrap_dir, pattern, excluded_tip)

  if(tree_style == "age"){
    if(time_file == ""){
      stop("The file containing the sampling time information is not provided!")
    }
    if(branch_num != 0){
      stop("The value in branch must be time in year!")
    }
    p = plot.tree.bootstrap.age(mytree, time_file, title, lextra, rextra, da)
  }else if(tree_style == "ci"){
    # write bootstrap tree to a nex file
    fbs = str_replace(tree_file, ".txt", ".bstrap.nex")
    # the tree has correct tip labels after using write.nexus
    write.nexus(mytree, file = fbs)

    if(branch_num == 0){
      if(time_file == ""){
        stop("The file containing the sampling time information is not provided!")
      }
      tree_ci = get.ci.tree(fbs, bstrap_dir2, labels, T, nex_pattern, excluded_tip = excluded_tip)
      p = plot.tree.ci.node(tree_ci, time_file, title, lextra, rextra, da, T, T)
    }else{
      if(scale_factor == 1){
        ci_prefix = "nmut_0.95_CI"
      }else{
        ci_prefix = "mutsize_0.95_CI"
      }
      tree_ci = get.ci.tree(fbs, bstrap_dir2, labels, T, nex_pattern, ci_prefix, excluded_tip = excluded_tip)

      if(ggtree_style == 0){
        p = plot.tree.ci.node.mut(tree_ci, title, lextra, rextra, da, T, T, scale_factor)
      }else{
        p = plot.tree.ci.node.mut.smpl(tree_ci, title, lextra, rextra, da, T, T, scale_factor)
      }
    }

    if(with_cn){
      d = fortify(tree_ci@phylo)
      ordered_nodes = d$label[order(d$y, decreasing = T)]
      d_seg = get.cn.data.by.pos(cn_file, pos_file, seg_file, cyto_file, labels, ordered_nodes, has_normal, bin_file, seed, is_haplotype_specific, cn_max, excluded_tip)
      # get the node order of the tree and reorder heatmap
      d_seg = d_seg %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))
      phmap = plot.cn.heatmap(d_seg, "")

      pc = ggarrange(p, phmap, nrow = 1, widths = c(6.5, 9.5))
      # pc = phmap %>% insert_left(p, width = 0.6)  # not work for internal nodes
      #  width = 16
      ggsave(out_file, pc, width = width, height = height)
    }else{
      ggsave(out_file, p, width = width, height = height)
    }

  }else{
    p = plot.tree.bootstrap(mytree, title, rextra, da, scale_factor)
    # width = 8
    ggsave(out_file, p, width = width, height = height)
  }

}else{
  message("plotting type not supported!")
}
