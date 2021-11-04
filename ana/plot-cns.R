#!/usr/bin/env Rscript

suppressMessages(library(copynumber))
suppressMessages(library(reshape))
suppressMessages(library(tools))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))


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
               panel.grid.minor=element_blank(),
               panel.grid.major=element_blank(),
               panel.spacing.y=unit(0, "lines"),
               panel.spacing.x=unit(0, "lines"),
               panel.border = element_rect(color = "#d8d8d8", fill = NA, size = .2),
               panel.background=element_rect(fill = "#f0f0f0")
)

cn_colors1 = c("#6283A9","#bdd7e7","#f0f0f0","#FCAE91", "#B9574E", "#76000D", "#8B0000", "#000000")
# Max CN to show in heatmap
max_cn = 7

get.cn.data <- function(d){
  md <- melt(d, id=c("sample","chromosome","index"))
  data <- cast(md, chromosome + index ~ sample)
  data %>% group_by(chromosome) %>% count() -> bin_count
  bins_4401 %>% group_by(chromosome) %>% count() -> ref_bin_count

  if(sum(bin_count$n) == sum(ref_bin_count$n)){
    # When the input data have 4401 bins
    data <- cbind(bins_4401, data)
  }
  else{ # When the data have less bins
    # Change reference bins to the desired number
    # adapted from https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html
    nested_bins <- bins_4401 %>%
      group_by(chromosome) %>%   # prep for work by Species
      nest() %>%              # --> one row per Species
      ungroup() %>%
      mutate(n=bin_count$n) # add sample sizes

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


# Assume in_file has absolute path
# Plot CNP separately
plot.cn <- function(in_file, out_file=""){
  dir <- dirname(in_file)
  stub <- file_path_sans_ext(basename(in_file), compression = T)
  if(out_file != ""){
    fout <- out_file
  }
  else{
    fout <- file.path(dir, paste("plot-",stub,".pdf",sep=""))
  }

  d <- read.table(in_file, header=FALSE)
  names(d) <- c("sample", "chromosome", "index", "cn")
  nsample <- length(unique(d$sample))

  data <- get.cn.data(d)

  #par(ask=F)
  pdf(fout, height=10, width=20)
  plotGenome(data=data, ylab="copy number", layout=c(2,3), sample=c(1:nsample), equalRange=FALSE, q=0, col="red" )
  dev.off()
}


# Plot CNPs of multiple samples as a heatmap for easy comparison
# Format of d_seg: sample, chrom, start, end, cn, ploidy
plot.cn.heatmap <- function(d_seg, fout, main, type="absolute", theme = theme1, cn_colors = cn_colors1, allele_specific = F, w = 12, h = 4){
    d_seg$cn = round(d_seg$cn)
    d_seg$chrom = factor(d_seg$chrom, levels=paste("",c(c(1:22)),sep=""))
    d_seg$pos = (d_seg$end + d_seg$start) / 2
    d_seg$width = (d_seg$pos - d_seg$start) * 2
    
    if(type=="absolute"){
      print("Plot absolute copy number")
      d_seg %>% mutate(cn=if_else(cn > max_cn, max_cn, cn)) -> d_seg
      if(!allele_specific){
        cn_vals = c("0", "1", "2", "3", "4", "5", "6", "7")
      }else{  # normal allele copy is 1
        cn_vals = c("-1", "0", "1", "2", "3", "4", "5", "6")
      }
    }else{
      print("Plot relative copy number")
      d_seg %>% mutate(cn = cn - ploidy) %>% mutate(cn=if_else(cn > max_cn, max_cn, cn)) %>% mutate(cn=if_else(cn < -2, -2, cn)) -> d_seg
      cn_vals = c("-2", "-1", "0", "1", "2", "3", "4", "5")
    }
    d_seg$cn = as.factor(d_seg$cn)
    # unique(d_seg$cn)
    
    plot<- ggplot(data=d_seg, aes(x=pos, y=sample, width=width)) +
      geom_tile(aes(fill=cn)) +
      facet_grid(sample ~ chrom, scales="free", space = "free", switch='y') +
      scale_color_manual(values=cn_colors, limits=cn_vals) +
      scale_fill_manual(values=cn_colors, limits=cn_vals) +
      scale_y_discrete(breaks=unique(d_seg$sample), expand = c(0,0))+
      scale_x_discrete(expand = c(0,0)) + theme
    plot = plot + ggtitle(main)
    
    ggsave(fout, plot, width = w, height = h)
}



option_list = list(
  make_option(c("-f", "--in_file"), type="character", default="",
              help="input file name [default=%default]", metavar="character"),
  make_option(c("-b", "--bin_file"), type="character", default="",
              help="bin file name [default=%default]", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default="",
              help="The name of output file [default=%default]", metavar="character"),
  make_option(c("-d", "--data_dir"), type="character", default="",
              help="The directory containing all the files to plot [default=%default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="",
              help="The naming pattern of files [default=%default]", metavar="character")
  # make_option(c("-t", "--plot_type"), type="character", default="single",
  #             help="The type of plot, including: all (plotting all files in a directory), single (plotting a single file) [default=%default]", metavar="character"),
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

in.file <- opt$in_file
data.dir <- opt$data_dir
bin.file <- opt$bin_file
out.file <- opt$out_file
pattern <- opt$pattern

plot_type <- "all"
if(in.file != ""){
  plot_type <- "single"
}


load(bin.file)

if(0){  # test
   d <- read.table("cnprofile-2.txt",header=FALSE)
   names(d) <- c("sample","chromosome","index","cn")
   nsample <- length(unique(d$sample))

   md <- melt(d, id=c("sample","chromosome","index"))
   data <- cast(md,chromosome+index~sample)
   data <- cbind(bins_4401, data)
   data <- data[,-c(3,4,5)]

   #par(ask=F)
   plotGenome(data=data, ylab="copy number", layout=c(2,3), sample=c(1:nsample), equalRange=FALSE,q=0, col="red" )
}

if (plot_type == "all"){
  dir <- data.dir
  cat(paste0("Plotting all related files in directory ", dir, "\n"))
  if(pattern == ""){
    files <- list.files(dir,"^sim\\-data\\-\\d+\\-cn")
  }
  else{
    files <- list.files(dir, pattern = pattern)
  }
  print(files)
  for(f in files){
    fname = file.path(dir, f)
    cat(paste0("Plotting file ", fname, "\n"))
    plot.cn(fname, out.file)
  }
}

if (plot_type == "single"){
  cat(paste0("Plotting the copy number variations in ", in.file, "\n"))
  plot.cn(in.file, out.file)
}
