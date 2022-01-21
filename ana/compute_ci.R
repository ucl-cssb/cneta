#!/usr/bin/env Rscript

# This script is used to compute confidence interval (CI) of a file with two columns where 1st column being ID and 2nd column being the data

args = commandArgs(trailingOnly = T)

fname = args[1]
digit = as.integer(args[2])
is_normal = args[3]

din = read.table(fname, strip.white = T, sep=" ")

data = din[,2]

# assume normal distribution
if(is_normal == "Y"){
  smu = mean(data)
  error = qnorm(0.975) * sd(data) / sqrt(length(data))
  left = smu - error
  right = smu + error
}else{
  # using 2.5th and 97.5th percentile
  sorted_data = sort(data)
  ci95 = quantile(sorted_data, probs = c(0.025, 0.975))
  left = ci95[1]
  right = ci95[2]
}


ci_str = paste0("(", round(left, digit), ", ", round(right, digit), ")")
cat(ci_str)
