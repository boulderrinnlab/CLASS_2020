######
# This script will run the full analysis using the singularity container 
#####

#### SETUP
library(rmarkdown)


# 01_consensus_peaks
rmarkdown::render("01_consensus_peaks/consensus_peak_set.Rmd", md_document(variant = "markdown_github"))

# 02_global_clustering
# - count_based_clustering.Rmd (bam files)
# - peak_based_clustering.Rmd (consensus peak files)
 
# 03_peak_feature_intersect
# - subset_peaks_by_feature_type.Rmd
# - permutation_test_by_feature_type.Rmd
 
# 04_correlation_with_expression
# - retrieve_expression.Rmd
# - tf_binding_vs_expression.Rmd
 
# 05_feature_clustering
# - feature_based_clustering.Rmd
# - clustering_per_feature_type.Rmd
# - promoter_profile.Rmd
 
# 06_repeat_analysis

