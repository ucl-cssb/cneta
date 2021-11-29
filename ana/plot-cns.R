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
