library(ape)
library(phangorn)
library(ggplot2)
library(tidyr)
# library(beeswarm)


if(!require(optparse)){
    install.packages("optparse")
    library(optparse)
}


##code from github sveta code
make.tree <- function(d){
   nedge <- nrow(d)
   nleaf <- (nedge + 2)/2
   nnode <- nleaf - 1

   mytree <- list()
   mytree$edge <- as.matrix(d[,c(1,2)])
   mytree$Nnode <- as.integer(nnode)
   mytree$tip.label <- paste(1:nleaf)
   mytree$edge.length <- d[,3]
   class(mytree) <- "phylo"
   checkValidPhylo(mytree)
   return(mytree)
}



compare.single.tree <- function(opt){
  ref_tree=opt$ref_tree
  inf_tree=opt$inf_tree
  file_plot=opt$file_plot
  file_dist=opt$file_dist

  # ref_tree="D:\\Gdrive\\git\\sveta\\example\\sim-data-1-tree.txt"
  # inf_tree="D:\\Gdrive\\git\\sveta\\example\\MaxL-m2-o1-sim-data-1.txt"
  # file_plot="D:\\Gdrive\\git\\sveta\\example\\blen-diff-1.pdf"
  # file_dist="D:\\Gdrive\\git\\sveta\\example\\tree-dist-1.txt"

  dd <- read.table(ref_tree, header=T)
  simtree <- make.tree(dd)

  dd <- read.table(inf_tree, header=T)
  inferred_tree <- make.tree(dd)

  # match the tree edges
  sim_edges = as.data.frame(simtree$edge)
  infer_edges = as.data.frame(inferred_tree$edge)
  sim_edges %>% unite(start, end, col = "edge",sep="_") -> sim_edges
  vedges_sim = as.vector(sim_edges$edge)
  infer_edges %>% unite(start, end, col = "edge",sep="_") -> infer_edges
  vedges_infer = as.vector(infer_edges$edge)
  idx_infer = match(vedges_sim, vedges_infer)

  total_diff = 0
  total_rel_diff = 0
  diff_blen = c()
  relative_diff_blen = c()
  for(j in 1:length(simtree$edge.length)){
    if(simtree$edge.length[j] == 0) next
    diff = inferred_tree$edge.length[idx_infer[j]]-simtree$edge.length[j]
    total_diff = total_diff + abs(diff)
    diff_blen = c(diff_blen, diff)

    rel_diff = diff / simtree$edge.length[j]
    relative_diff_blen = c(relative_diff_blen, rel_diff)
    total_rel_diff = total_rel_diff + abs(rel_diff)
  }
  # diff_blen

  avg_diff = total_diff/length(simtree$edge.length)
  avg_rel_diff = total_rel_diff/length(simtree$edge.length)

  # par(mfrow=c(2,1))
  # plot.tree(simtree)
  # plot.tree(inferred_tree)

  pdf(file_plot)
  # Print bar plot of the estimation errors
  barplot(diff_blen)
  barplot(relative_diff_blen)
  dev.off()


  tree_dist_all=file.path(file_dist)
  close(file( tree_dist_all, open="w" ))
  header = paste("total_blen_diff", "avg_blen_diff", "total_relative_diff", "avg_relative_diff", "normalizedRF", "RF", "branch", "path", "weightedPath", sep="\t")
  write(header, tree_dist_all, append=F)

  tdist <- treedist(simtree, inferred_tree)
  nRFdist <- RF.dist(simtree, inferred_tree, normalize=T)
  line=paste(total_diff, avg_diff, total_rel_diff, avg_rel_diff, nRFdist, tdist[1], tdist[2], tdist[3], tdist[4], sep="\t")
  write(line, tree_dist_all, append=T)
}


compare.all.trees <- function(opt){
    # res_all = data.frame(id = integer(), RF_dist = numeric(), blen_dist = numeric(), num_mut = integer(), num_mut_tip = integer(), mut_diff = numeric());
    res_all = data.frame()

    fml=file.path(opt$tree_dir, "ml_finished")
    lnl_all = read.table(fml)
    names(lnl_all) = c("id", "neg_lnl")
    head(lnl_all)
    str(lnl_all)

    sample_dlarge = c()
    # Trees for which real tree and best tree do not match
    fmis=file.path(opt$tree_dir, "tree_mismatch")
    if(file.exists(fmis) && file.info(fmis)$size > 0){
      tree_mismatch = read.table(fmis)
      tree_mismatch = tree_mismatch$V1
    }else{
      tree_mismatch = NULL
    }


    tree_dlarge=file.path(opt$tree_dir, "tree_dist_large")
    close( file( tree_dlarge, open="w" ) )

    tree_dist_all=file.path(opt$tree_dir, "tree_dist_all")
    close( file( tree_dist_all, open="w" ) )
    header = paste("id", "normalizedRF", "RF", "branch", "path", "weightedPath", sep="\t")
    write(header, tree_dist_all, append=T)


    pdf(opt$file_plot)

    # maxLtreefiles <- list.files(pattern="plot-MaxL")
    # maxLtreefiles
    for (i in 1:opt$nsim){
    # for (k in 1:length(maxLtreefiles)){
        # file = maxLtreefiles[k]
        # fname = gsub(pattern = "\\.pdf$", "", file)
        # fields = strsplit(fname, "-")
        # i = fields[[1]][length(fields[[1]]) - 1]
        # j = fields[[1]][length(fields[[1]])]
        # i = 1
        if(opt$incl_all==0){
          if(i %in% tree_mismatch) next
        }

        simtreename <- file.path(opt$data_dir, paste0("sim-data-", i, "-tree.txt"))
        maxLtreename1 <- file.path(opt$tree_dir, paste0("MaxL-m", opt$model, "-o1-sim-data-", i, ".txt"))
        maxLtreename2 <- file.path(opt$tree_dir, paste0("MaxL-o1-sim-data-", i, ".txt"))

        if(file.exists(maxLtreename1)){
          maxLtreename = maxLtreename1
        }else{
          maxLtreename = maxLtreename2
        }

        # Some results may not exist due to failure in WGD introduction
        if(!file.exists(maxLtreename)) next

        run <- paste0("run", i)
        dd <- read.table(simtreename, header=T)
        simtree <- make.tree(dd)

        num_mut <- sum(dd$nmut)
        num_tip = length(simtree$tip.label)
        tip <- dd[dd$end < num_tip,]
        num_mut_tip <- sum(tip$nmut)
        num_sample=num_tip - 1

        tip_comb = combn(tip$end, 2)
        sum = 0
        for(j in 1:ncol(tip_comb)){
            t1 = tip_comb[1,j]
            t2 = tip_comb[2,j]
            m1 = tip[tip$end==t1,]$nmut
            m2 = tip[tip$end==t2,]$nmut
            diff = abs(m1-m2)
            sum = sum + diff
        }
        mut_diff = sum / ncol(tip_comb)
        dd <- read.table(maxLtreename, header=T)
        maxLtree <- make.tree(dd)

        tdist <- treedist(simtree, maxLtree)
        nRFdist <- RF.dist(simtree, maxLtree, normalize=T)

        if(tdist[1] > 0){
          write(i, tree_dlarge, append=T)
          print(paste0("Distance for ", i , " is ", tdist[1]))
          sample_dlarge = c(sample_dlarge, i)
        }
        res = data.frame(id = i, RF_dist=as.numeric(tdist[1]), nRF_dist=as.numeric(nRFdist), blen_dist=as.numeric(tdist[2]), num_sample=num_sample, num_mut=num_mut, num_mut_tip = num_mut_tip, mut_diff = mut_diff)
        res_all <- rbind(res_all, res)

        # Compare the estimation of branch lengths (differences between real value and inferred value)
        diff_blen = c()
        for(j in 1:length(simtree$edge.length)){
          diff = simtree$edge.length[j] - maxLtree$edge.length[j]
          diff_blen = c(diff_blen, diff)
        }

        # Print bar plot of the estimation errors
        barplot(diff_blen)

        line=paste(i, nRFdist, tdist[1], tdist[2], tdist[3], tdist[4], sep="\t")
        write(line, tree_dist_all, append=T)
    }

    head(res_all)
    str(res_all)


    all = merge(res_all, lnl_all, by=c("id"))
    head(all)
    str(all)



    # pdf("dist_hist.pdf")
    p1 = ggplot(res_all, aes(x=RF_dist)) + geom_histogram(binwidth=1) + stat_bin(binwidth=1, geom="text", size=3.5, aes(label=..count..), vjust=-1.5)
    # dev.off()
    print(p1)


    p1 = ggplot(res_all, aes(x=num_sample, y=nRF_dist)) + geom_violin() + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=23, size=2, color="red")

    # dev.off()
    print(p1)


    # pdf("distVmutdiff.pdf")
    p2 = ggplot(res_all, aes(x=as.factor(RF_dist), y=mut_diff)) + geom_boxplot() + geom_jitter()
    # dev.off()
    print(p2)

    # pdf("distVlnl.pdf")
    p3 = ggplot(all, aes(x=as.factor(RF_dist), y=neg_lnl)) + geom_boxplot() + geom_jitter()
    # dev.off()
    print(p3)

    p4 = ggplot(res_all, aes(x=as.factor(RF_dist), y=num_mut)) + geom_boxplot() + geom_jitter()
    print(p4)

    dev.off()
    #write(MutDist, "MutDist.txt")

    # Error <- matrix(0, 2, 2, dimnames = list(c("Run", "Tree building"), c("Number", "%")))
    # Error[1,2] <- length(maxLtreefiles)/length(treefiles)
    # Error[1,1] <- length(maxLtreefiles)
    # RFDistFreq <- table(distance[,1])
    # Error[2,2] <- RFDistFreq[1]/sum(RFDistFreq)
    # Error[2,1] <- RFDistFreq[1]
    #
    # write.table(distance, "NewMutDist.txt", row.names=TRUE)
    # write.table(Error, "Error.txt", col.names=T)
    # write.table(RFDistFreq, "DistRFFreq.txt", row.names=FALSE)

    ##add in the boxplots, but the distance file has all the information that is needed
    ##also add in error rate of distance and error rate of re-constructung the trees

}


option_list = list(
make_option(c("-t", "--tree_dir"), type="character", default="",
              help="The directory containing all the files of tree building programs [default=%default]", metavar="character"),
make_option(c("-d", "--data_dir"), type="character", default="",
              help="The directory containing all the simulated files [default=%default]", metavar="character"),
make_option(c("-r", "--ref_tree"), type="character", default="",
              help="The reference tree [default=%default]", metavar="character"),
make_option(c("-i", "--inf_tree"), type="character", default="",
              help="The inferred tree [default=%default]", metavar="character"),
make_option(c("-o", "--file_plot"), type="character", default="stat_plot.pdf",
              help="The name of output plot file [default=%default]", metavar="character"),
make_option(c("-s", "--file_dist"), type="character", default="tree_dist.txt",
            help="The name of output file [default=%default]", metavar="character"),
make_option(c("-a", "--incl_all"), type="integer", default=0,
              help="Whether or not to exclude tree for which whose maximum likelihood is smaller than the maximum likelihood of all trees [default=%default]", metavar="integer"),
make_option(c("-m", "--model"), type="integer", default=2,
            help="The model of evolution [default=%default]", metavar="integer"),
make_option(c("-p", "--type"), type="integer", default=0,
            help="The type of analysis (0: multiple tree files in a directory; 1: a single tree) [default=%default]", metavar="integer"),
make_option(c("-n", "--nsim"), type="integer", default=100,
              help="The number of simulations [default=%default]", metavar="integer")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(opt$type == 0){
  compare.all.trees(opt)
}else{
  compare.single.tree(opt)
}
