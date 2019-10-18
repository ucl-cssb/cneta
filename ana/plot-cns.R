library(copynumber)
library(reshape)
library(tools)
library(tidyr)
library(dplyr)
library(purrr)

if(!require(optparse)){
    install.packages("optparse")
    library(optparse)
}

# load("../../lpWGS phylogenies/bin_locations_4401.Rdata")
#load("./bin_locations_4401.Rdata")

# args = commandArgs(trailingOnly=TRUE)
# data.dir <- args[1]
# bin.file <- args[2]
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
if(in.file!=""){
  plot_type <- "single"
}

load(bin.file)

# Assume in_file has absolute path
print.cn <- function(in_file, out_file=""){
  dir <- dirname(in_file)
  stub <- file_path_sans_ext(basename(in_file), compression = T)
  if(out_file!=""){
    fout <- out_file
  }
  else{
    fout <- file.path(dir, paste("plot-",stub,".pdf",sep=""))
  }

  d <- read.table(in_file, header=FALSE)
  names(d) <- c("sample","chromosome","index","cn")
  nsample <- length(unique(d$sample))

  md <- melt(d, id=c("sample","chromosome","index"))
  data <- cast(md,chromosome+index~sample)
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

  #par(ask=F)
  pdf(fout, height=10,width=20)
  plotGenome(data=data, ylab="copy number", layout=c(2,3), sample=c(1:nsample), equalRange=FALSE,q=0, col="red" )
  dev.off()
}


if(0){
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
  #dir <- "../sim-data/"
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
    print.cn(fname, out.file)
  }
}

if (plot_type == "single"){
  cat(paste0("Plotting the copy number variations in ", in.file, "\n"))
  print.cn(in.file, out.file)
}
