#!/usr/bin/env Rscript

# This script is used to convert the simulated copy numbers to the format that can be accepted by MEDICC2
# Sample command: Rscript convert4medicc2.R -f sim-data-1-allele-cn.txt.gz -m medicc_input.tsv

library(readr)
library(dplyr)
library(optparse)


option_list = list(
  make_option(c("-i", "--file_cn"), type = "character", default = "",
              help = "input haplotype-specific copy number file [default=%default]", metavar = "character"),
  make_option(c("-o", "--file_out"), type = "character", default = "",
              help = "The output TSV file to be used as input to MEDICC2 [default=%default]", metavar = "character")
);

opt = parse_args(OptionParser(option_list = option_list))


fin = opt$file_cn
fout = opt$file_out

ac = read.table(fin)

names(ac) = c("sid", "chr", "seg", "cnA", "cnB")

# exclude normal sample
normal_id = max(unique(ac$sid))

# Find major and minor copy
ac %>% filter(sid != normal_id) %>% rowwise() %>% mutate(major = max(cnA, cnB), minor = min(cnA, cnB)) -> ac_full

ac_full %>% select(sample_id = sid, chrom = chr, start = seg, end = seg, cn_a = major, cn_b = minor) %>% mutate(end = end + 1, sample_id = paste0("taxon_", sample_id), chrom = paste0("chrom", chrom)) -> cn_medicc

write_tsv(cn_medicc, fout)