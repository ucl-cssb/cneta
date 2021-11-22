#!/usr/bin/env Rscript

suppressMessages(library(copynumber))
suppressMessages(library(reshape))
suppressMessages(library(tools))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))


get.cn.data.by.bin <- function(d, bins_4401){
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



get.outfile.name <- function(in_file){
  dir <- dirname(in_file)
  stub <- file_path_sans_ext(basename(in_file), compression = T)

  fout <- file.path(dir, paste("plot-",stub,".pdf",sep=""))

  return(fout)
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



option_list = list(
  make_option(c("-f", "--in_file"), type="character", default="",
              help="The input file which contains the copy numbers [default=%default]", metavar="character"),
  make_option(c("-b", "--bin_file"), type="character", default="",
              help="The file which contains the location of each bin when the bin size is 500 Kbp [default=%default]", metavar="character"),
  # make_option(c("", "--pos_file"), type="character", default="",
              # help="The file which contains the positions of each site [default=%default]", metavar="character"),  
  make_option(c("-o", "--out_file"), type="character", default="",
              help="The name of output file [default=%default]", metavar="character"),
  make_option(c("-d", "--data_dir"), type="character", default="",
              help="The directory containing all the files to plot [default=%default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="",
              help="The naming pattern of files [default=%default]", metavar="character")
  # make_option(c("", "--plot_style"), type="character", default="line",
               # help="The style of plot, including: line, heatmap [default=%default]", metavar="character"),
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

in_file <- opt$in_file
data_dir <- opt$data_dir
bin_file <- opt$bin_file
out_file <- opt$out_file
pattern <- opt$pattern

plot_type <- "all"
if(in_file != ""){
  plot_type <- "single"
}

load(bin_file)


if(plot_type == "all"){
  dir <- data_dir
  cat(paste0("Plotting all related files in directory ", dir, "\n"))
  if(pattern == ""){
    files <- list.files(dir,"^sim\\-data\\-\\d+\\-cn")
  }else{
    files <- list.files(dir, pattern = pattern)
  }
  print(files)
  for(f in files){
    fname = file.path(dir, f)
    cat(paste0("Plotting file ", fname, "\n"))
    out_file = get.outfile.name(fname)
    plot.cn(fname, out_file, bins_4401)
  }
}else{    # plot_type == "single"
  cat(paste0("Plotting the copy number variations in ", in_file, "\n"))
  if(out_file == ""){
    out_file = get.outfile.name(in_file)
  }
  plot.cn(in_file, out_file, bins_4401)
}
