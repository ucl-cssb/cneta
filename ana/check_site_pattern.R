library(ggplot2)
library(tidyr)
library(dplyr)
library(optparse)


# Check site patterns in the input copy number data
# Output: the count of each unique pattern across all samples

option_list = list(
make_option(c("-c", "--file_cn"), type="character", default="",
              help="The input file containing the copy numbers with four colomns ("sample", "chr", "seg", "tcn") [default=%default]", metavar="character"),
make_option(c("-t", "--file_txt"), type="character", default="site_pattern.txt",
              help="The name of output file (in txt format) [default=%default]", metavar="character")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);


in_file=opt$file_cn
d <- read.table(in_file, header=FALSE)
names(d)=c("sample", "chr", "seg", "tcn")

d %>% unite(site, c(chr, seg)) -> du
length(unique(du$site))
du %>% spread(site, tcn) -> dp

tdp=t(dp)
colnames(tdp)=tdp[1,]
tdp = as.data.frame(tdp[2:nrow(tdp),])
tdp %>% group_by_all() %>% summarise(COUNT = n())-> sp
write.table(sp, opt$file_txt, row.names = F, col.names = T, sep="\t")
