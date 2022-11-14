#!/usr/bin/env Rscript

library(tidyverse)
library(MCMCpack)

# This script is used to create data with subclonal structures based on the output of CNETS (only applicable to total copy number data)


# assume the last sample is normal, R pass function by value
# simulate data with normal mixture
# purs: purities of each sample
convert_cns_by_purity <- function(cns_sel, purs){
  for(i in 1:nrow(cns_sel)){
    cns_sel[i,] = round(purs[i] * 2 + (1 - purs[i]) * cns_sel[i,])
  }
  return(cns_sel)
}


# get clone IDs in a vector
get_clones_with_primary <- function(i, nclone, nprop, incl_normal, idx_normal_clone){
  clones = c()
  if(incl_normal){
    clones = sample(1:nclone, nprop - 2)
    while(i %in% clones){
      clones = sample(1:nclone, nprop - 2)
    }
    clones = c(i, clones, idx_normal_clone)
  }else{
    clones = sample(1:nclone, nprop - 1)
    clones = c(i, clones)
  }
  
  return(clones)
}


# change the copy numbers of each clone/sample by mixing with other clones/samples
# assume several clones with different proportions, with the proportion of the first clone fixed and the last clone being normal
convert_cns_by_props <- function(cns_sel, nprop = 3, incl_normal = T, orig_props = c(), major_prop = 0){
  nclone = nrow(cns_sel)
  cns_new = matrix(0, nrow = nrow(cns_sel), ncol = ncol(cns_sel))
  colnames(cns_new) = colnames(cns_sel)
  idx_normal_clone = nrow(cns_sel)
  if(incl_normal){
    nclone = nclone - 1
    cns_new[idx_normal_clone, ] = 2
  }
  df_prop_clone = data.frame()
  
  # assume each row is a sample in the new matrix 
  # assume each row is a clone in the original matrix
  # randomly pick up nprop clones 
  # nclone is the number of tumour clones
  for(i in 1:nclone){
    if(major_prop > 0){
      alpha = rep(1, nprop - 1)
      props = t(rdirichlet(1, alpha))
      props = props * (1 - major_prop)
      props = c(major_prop, props)
      # ensure the last clone is normal
      clones = get_clones_with_primary(i, nclone, nprop, incl_normal, idx_normal_clone)
       
    }else{
      if(length(orig_props) == 0){
        alpha = rep(1, nprop)
        props = t(rdirichlet(1, alpha))
        clones = c()
        # ensure the last clone is normal
        if(incl_normal){
          clones = sample(1:nclone, nprop - 1)
          clones = c(clones, idx_normal_clone)
        }else{
          clones = sample(1:nclone, nprop)
        }
      }else{
        props = orig_props
        # ensure first clone is the same
        clones = get_clones_with_primary(i, nclone, nprop, incl_normal, idx_normal_clone)
      }     
    }

    # print(props)
    
    df = data.frame(sample=i)
    df_clone = as.data.frame(t(clones))   # weird to need t()
    names(df_clone) = paste0("clone", 1:nprop)
    df_prop = as.data.frame(t(props))
    names(df_prop) = paste0("prop", 1:nprop)
    df = cbind(df, df_clone, df_prop)
    df_prop_clone = rbind(df_prop_clone, df)
    
    for(j in 1:nprop){
      cns_new[i,] = cns_new[i,] + cns_sel[clones[j],] * props[j]
    }
    cns_new[i,] = round(cns_new[i,])
  }
  
  return(list(df_prop_clone = df_prop_clone, cns_new = cns_new))
}



test_conversion <- function(cns_mat, purs){
  # test on a small matrix
  cols = sample(1000, 10)
  cns_sel = cns_mat[, cols]
  cns_orig = cns_sel
  cns_sel = convert_cns_by_purity(cns_sel, purs)
  cns_sel
  cns_orig
  cns_sel - cns_orig
}


test_conversion_by_props <- function(cns_mat, nprop){
  # test on a small matrix
  cols = sample(1000, 10)
  cns_sel = cns_mat[, cols]
  cns_orig = cns_sel
  cns_sel = convert_cns_by_props(cns_sel, nprop)
  cns_sel
  cns_orig
  cns_sel - cns_orig
}


# file_cn: the total copy number file simulated by CNETS
# the purity of each sample is drawn from a normal distribution
convert_one_file_purity <- function(file_cn, odir, mu_purity = 0.5, sd_purity = 0.3){
  data_cn <- read.table(file_cn)
  names(data_cn) <- c("sid", "chr", "seg", "cn")
  
  file_base = basename(file_cn)
  
  # transformation
  data_cn %>%
    unite(chr, seg, col = "pos", sep = "_") %>%
    spread(pos, cn) %>%
    dplyr::select(-sid)  -> cns_wide
  cns_mat = as.matrix(cns_wide)
  
  # convert with purity alone
  purity_mean = mu_purity
  purity_sd = sd_purity
  nclone = nrow(cns_wide) - 1
  purs = rnorm(nclone, purity_mean, purity_sd)
  purs[purs > 1] = 1
  
  fname = str_replace(file_base, "-cn.txt.gz", "-purity.txt")
  fout = file.path(odir, fname)
  # fout = file.path(odir, paste0("sim-data-", i, "-purity.txt"))
  write_tsv(as.data.frame(purs), fout, col_names = F)
  
  purs = c(purs, 1)
  
  cns_with_purity = convert_cns_by_purity(cns_mat, purs)
  # diff = cns_mat - cns_with_purity
  # which(diff != 0)
  
  # final output: same format as the input
  pos_names = colnames(cns_with_purity)
  cns_with_purity = cns_with_purity %>% as.data.frame()
  cns_with_purity$sid <- 1:nrow(cns_with_purity)
  cns_with_purity %>% gather(all_of(pos_names), key = "pos", value = "cn") %>%
    separate(pos, sep = "_", into = c("chr", "seg")) %>% rowwise() %>% 
    mutate(chr = as.integer(chr), seg = as.integer(seg)) %>%
    arrange(sid, chr, seg) -> dcns_new
  
  # write the data into gz file for tree building
  # fout = file.path(odir, paste0("sim-data-", i, "-cn.txt.gz"))
  fout = file.path(odir, file_base)
  gz1 <- gzfile(fout, "w")
  write.table(dcns_new, gz1, quote = F, row.names = F, col.names = F, sep = "\t")
  close(gz1)
}


# file_cn: the total copy number file simulated by CNETS
# nprop: number of clones in each sample (including normal clone)
convert_one_file_subclone <- function(file_cn, odir, nprop = 3, incl_normal = T, props = c(), major_prop = 0){
  data_cn <- read.table(file_cn)
  names(data_cn) <- c("sid", "chr", "seg", "cn")
  
  file_base = basename(file_cn)
  
  # transformation
  data_cn %>%
    unite(chr, seg, col = "pos", sep = "_") %>%
    spread(pos, cn) %>%
    dplyr::select(-sid)  -> cns_wide
  
  cns_mat = as.matrix(cns_wide)

  res = convert_cns_by_props(cns_mat, nprop, incl_normal, props, major_prop)
  
  df_prop_clone = res$df_prop_clone
  cns_with_subclone = res$cns_new
  # diff = cns_mat - cns_with_subclone
  # which(diff != 0)
  
  # final output: same format as the input
  pos_names = colnames(cns_with_subclone)
  cns_with_subclone = cns_with_subclone %>% as.data.frame()
  cns_with_subclone$sid <- 1:nrow(cns_with_subclone)
  cns_with_subclone %>% tidyr::gather(all_of(pos_names), key = "pos", value = "cn") %>%
    tidyr::separate(pos, sep = "_", into = c("chr", "seg")) %>% dplyr::rowwise() %>% 
    dplyr::mutate(chr = as.integer(chr), seg = as.integer(seg)) %>%
    dplyr::arrange(sid, chr, seg) -> dcns_new
  
  # write the data into gz file for tree building
  # fout = file.path(odir, paste0("sim-data-", i, "-cn.txt.gz"))
  fout = file.path(odir, file_base)
  gz1 <- gzfile(fout, "w")
  write.table(dcns_new, gz1, quote = F, row.names = F, col.names = F, sep = "\t")
  close(gz1)
  
  # fout = file.path(odir, paste0("sim-data-", i, "-prop.txt"))
  fname = str_replace(file_base, "-cn.txt.gz", "-prop.txt")
  fout = file.path(odir, fname)
  write_tsv(df_prop_clone, fout, col_names = F)
}






