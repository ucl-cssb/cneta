#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))


# Check site patterns in the input copy number data
# Output: the count of each unique pattern across all samples

option_list = list(
make_option(c("-c", "--file_cn"), type="character", default="",
              help="The input file containing the copy numbers with four colomns (\"sample\", \"chr\", \"seg\", \"tcn\") [default=%default]", metavar="character"),
make_option(c("-t", "--file_txt"), type="character", default="site_pattern.txt",
              help="The name of output file in txt format; row: site pattern (CN for each sample) followed by the number of this pattern [default=%default]", metavar="character")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

in_file = opt$file_cn
out_file = opt$file_txt

d <- read.table(in_file, header=FALSE)
names(d)=c("sample", "chr", "seg", "tcn")

d %>% unite(site, c(chr, seg)) -> du
du %>% spread(site, tcn) -> dp

tdp = t(dp)
colnames(tdp) = tdp[1,]
tdp = as.data.frame(tdp[2:nrow(tdp),])
group_by_all(tdp) %>% tally() -> sp

write.table(sp, out_file, row.names = F, col.names = T, sep = "\t")

# get maximum state changes
last = max(d$sample)
sp$max_site_change <- NA
ploidy = 2
for (i in 1:NROW(sp)) {
  max_change <- max(abs(sp[i, 1:last] - ploidy))
  sp$max_site_change[i] <- max_change
}
# print max site change to standard output for pipeline
cat(max(sp$max_site_change))
