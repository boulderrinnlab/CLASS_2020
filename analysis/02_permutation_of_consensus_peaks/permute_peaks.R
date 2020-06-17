knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(GenomicRanges)
library(rslurm)
library(regioneR)
library(effectsize)
source("../util/intersect_functions.R")
source("../util/_setup.R")

# IMPORT peaks and features to test overlaps during the randomizations(permutations) ####
peak_list <- import_peaks("../01_consensus_peaks/results/consensus_peaks/filtered_by_peaks/")

lncrna_promoters <- rtracklayer::import("../01_consensus_peaks/results/lncrna_promoters.gtf")
mrna_promoters <- rtracklayer::import("../01_consensus_peaks/results/mrna_promoters.gtf")

lncrna_matrix <- count_peaks_per_feature(lncrna_promoters, peak_list, type = "occurrence")
mrna_matrix <- count_peaks_per_feature(mrna_promoters, peak_list, type = "occurrence")

rmsk <- import_repeatmasker()
rmsk_family <- subset_rmsk(rmsk, rep_level = "family")
names(rmsk_family) <- paste(names(rmsk_family), "family", sep = "_")
rmsk_class <- subset_rmsk(rmsk, rep_level = "class")
names(rmsk_class) <- paste(names(rmsk_class), "class", sep = "_")

hg38 <- getGenome("hg38")

region_list <- c("lncrna_promoters" = list(lncrna_promoters), 
                 "mrna_promoters" = list(mrna_promoters), 
                 rmsk_family, rmsk_class)


canonical_chr <- as.character(unique(seqnames(region_list[[1]])))
# Sanitize RMSK to only those repeats on canonical chromosomes
for(i in 1:length(region_list)) {
  region_list[[i]] <- region_list[[i]][which(seqnames(region_list[[i]]) %in% canonical_chr)]
}

#### DEFINE Permutation TESTS ####
pars <- expand.grid(1:length(region_list), 1:length(peak_list)) %>% 
  as.data.frame()
names(pars) <- c("region_index", "peak_index")



#### TEST FUNCTION ####
perm_test <- function(region_index, peak_index, npermutations = 1000) {
  
  set.seed(12044593)
  region <- names(region_list)[[region_index]]
  tf <- names(peak_list)[[peak_index]]
  
  cat("Running overlap test for ", region, "  & ", tf, "\n\n")
  a_regions <- region_list[[region_index]]
  b_regions <- peak_list[[peak_index]]
  
  suppressWarnings(pt <- overlapPermTest(A = a_regions, 
                                         B = b_regions, 
                                         ntimes = npermutations, 
                                         non.overlapping = FALSE, 
                                         verbose = FALSE,
                                         genome = hg38,
                                         alternative =  "auto", 
                                         mc.cores = 1))
  
  ptdf <- data.frame("region" = region,
                     "tf" = tf,
                     "pval" = pt$numOverlaps$pval,
                     "zscore" = pt$numOverlaps$zscore,
                     "nperm" = pt$numOverlaps$ntimes,
                     "alternative" = pt$numOverlaps$alternative,
                     "observed" = pt$numOverlaps$observed,
                     "permuted" = paste(pt$numOverlaps$permuted, collapse = ";"))
  return(ptdf)
}

# Note this took about 400 hours to run
# TODO: re-write to be faster
res_files <- list.files("_rslurm_perm_overlaps/", full.names = T, pattern = "results")
if(length(res_files) == 0) {
  sjob <- slurm_apply(perm_test, pars, jobname = 'perm_overlaps',
                      add_objects = c("region_list", "peak_list", "hg38", "overlapPermTest"),
                      nodes = 22, cpus_per_node = 30, 
                      slurm_options = list(time = '400:00:00'),
                      submit = FALSE)
}
