In this folder, there are some R scripts and Python scripts for plotting, preprocessing, and testing.

bin_locations_4401.Rdata contains the position of 4401 bins in the human reference genome hg19.
The bins are obtained by non-overlapping sliding windows of size 500,000 bp along each chromosome.

For plotting of trees with copy number heat map, 

get.cn.data.by.pos(cn_file, pos_file, seg_file, cyto_file, labels, ordered_nodes, has_normal, bin_file, seed, is_haplotype_specific)

cn_file: 5 columns, with cnA, cnB
pos_file: input copy number, 5 columns
seg_file: 1 value for each sample, state, need to change to total CN