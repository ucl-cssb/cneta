#!/usr/bin/env Rscript

# This script is used to compute confidence interval (CI) of a file with two columns where 1st column being ID and 2nd column being the data

args = commandArgs(trailingOnly = T)

fname = args[1]
digit = as.integer(args[2])

din = read.table(fname, strip.white = T, sep=" ")

data = din[,2]

smu = mean(data)
error = qnorm(0.975) * sd(data) / sqrt(length(data))
left = smu - error
right = smu + error

ci_str = paste0("(", round(left, digit), ", ", round(right, digit), ")")
cat(ci_str)
