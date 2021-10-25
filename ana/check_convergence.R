#!/usr/bin/env Rscript

# Check the convergence of MCMC chains
# https://github.com/danlwarren/RWTY

library(rwty)

args = commandArgs(trailingOnly=TRUE)
dir = args[1]
fout = args[2]

# dir = "D:/Gdrive/git/sveta/example/mcmc"
# fout = "D:/Gdrive/git/sveta/example/mcmc/mcmc_rwty.pdf"



# tree_file = file.path(dir, "mcmc-sim-data-1_2.t")
# log_file = file.path(dir, "mcmc-sim-data-1_2.p")
#
# # Load one chain
# # chain1 <- load.trees(tree_file, type = "nexus", logfile = log_file, skip = 0, gens.per.tree = 1)
# chain1 <- load.trees(tree_file, logfile = log_file, format = "mb")
# str(chain1)
# chain1.rwty <- analyze.rwty(chain1, burnin=0, fill.color = 'lnl')
# # to see which plots you have
# names(chain1.rwty)
# makeplot.all.params(chain1.rwty, burnin=0)
# approx.ess <- topological.approx.ess(chain1, burnin = 0)
#


# Load multiple chains
multi.trees <- load.multi(dir, format = "mb")
str(multi.trees)
multi.rwty <- analyze.rwty(multi.trees, burnin=0, fill.color = 'lnl', filename = fout, overwrite = T)
names(multi.rwty)
# print(multi.rwty)

# makeplot.all.params(multi.rwty, burnin=0)
# approx.ess <- topological.approx.ess(multi.trees, burnin = 0)
# multi.rwty$lnl.trace$trace.plot
# multi.rwty$lnl.trace$density.plot
