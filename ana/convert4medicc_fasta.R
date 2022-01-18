#!/usr/bin/env Rscript

# This script is used to convert the simulated copy numbers to the format that can be accepted by MEDICC (copied from Simone)
# Sample command: Rscript convert4medicc.R -f sim-data-1-allele-cn.txt -d medicc_input

library(tidyr)
library(dplyr)
library(seqinr)
library(optparse)


option_list = list(
  make_option(c("-i", "--file_cn"), type = "character", default = "",
              help = "input haplotype-specific copy number file [default=%default]", metavar = "character"),
  make_option(c("-d", "--dir_out"), type = "character", default = "",
              help = "The directory to store results [default=%default]", metavar = "character")
);

opt = parse_args(OptionParser(option_list = option_list))

dir = opt$dir_out
ac_name = opt$file_cn
ac = read.table(ac_name)

names(ac) = c("sid", "chr", "seg", "cnA", "cnB")

dir.create(dir, showWarnings = F)

# Find major and minor copy
ac %>% rowwise() %>% mutate(major = max(cnA, cnB), minor = min(cnA, cnB)) -> ac_full
# t=ac_full$major-ac_full$cnA
# ac_full[which(t>0),]
# t=ac_full$minor-ac_full$cnA
# ac_full[which(t<0),]


chrs = unique(ac$chr)
# Write description file
file.desc = file.path(dir, paste0("desc.txt"))
close(file(file.desc, open = "w" ))

sample_prefix = "taxon_"
ncol = 1e6

for (i in 1:length(chrs)) {
  c = chrs[i]
  print(c)
  ac_full_chr = filter(ac_full, chr == c)
  # Find number of segments for each sample
  nsegs = table(ac_full_chr$sid)
  suffix = c

  file.out = file.path(dir, paste0("major_chr", suffix, ".fasta"))
  cn_major = select(ac_full_chr, c("sid", "major"))
  cn_major$id = rep(seq.int(nsegs[1]), dim(nsegs))
  cn_major_table = spread(cn_major, key = "id", value = "major")
  cn_major_seq = list()
  for (j in nrow(cn_major_table):1) {
    if (j == nrow(cn_major_table)) {
      sname = "diploid"
    }
    else {
      sname = paste0(sample_prefix, cn_major_table[j,]$sid)
    }
    cn_major_seq[[sname]] = as.character(cn_major_table[j, -1])
  }
  write.fasta(sequences = cn_major_seq, names = names(cn_major_seq), nbchar = ncol, file.out = file.out)


  file.out = file.path(dir, paste0("minor_chr", suffix, ".fasta"))
  cn_minor = select(ac_full_chr, c("sid", "minor"))
  cn_minor$id = rep(seq.int(nsegs[1]), dim(nsegs))
  cn_minor_table = spread(cn_minor, key = "id", value = "minor")
  cn_minor_seq = list()
  for (j in nrow(cn_minor_table):1) {
    if (j == nrow(cn_minor_table)) {
      sname = "diploid"
    }
    else{
      sname = paste0(sample_prefix, cn_minor_table[j,]$sid)
    }
    cn_minor_seq[[sname]] = as.character(cn_minor_table[j, -1])
  }
  write.fasta(sequences = cn_minor_seq, names = names(cn_minor_seq), nbchar = ncol, file.out = file.out)

  line = paste(paste0("chrom", c),  paste0("major_chr", suffix, ".fasta"), paste0("minor_chr", suffix, ".fasta"))
  write(line,file = file.desc,append = TRUE)
}
