######
# This script will run the full analysis using the singularity container 
#####

#### SETUP
library(rmarkdown)


# 01_consensus_peaks
rmarkdown::render("01_consensus_peaks/01_consensus_peaks.Rmd", md_document(variant = "markdown_github"))


# 02_permutation_of_consensus_peaks
rmarkdown::render("02_permutation_of_consensus_peaks/permutation_feature_intersects.Rmd", md_document(variant = "markdown_github"))
  

# 03_global_clustering
rmarkdown::render("03_global_clustering/03_global_clustering.Rmd", md_document(variant = "markdown_github"))
