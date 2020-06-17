######
# This script will run the full analysis using the singularity container 
#####

#### SETUP
library(rmarkdown)


# 01_consensus_peaks
rmarkdown::render(file.path("01_consensus_peaks",
                            "01_consensus_peaks.Rmd"), 
                  md_document(variant = "markdown_github"))


# 02_permutation_of_consensus_peaks
rmarkdown::render(file.path("02_permutation_of_consensus_peaks",
                            "permutation_feature_intersects.Rmd"), 
                  md_document(variant = "markdown_github"))
  

# 03_global_clustering
rmarkdown::render(file.path("03_global_clustering",
                            "03_global_clustering.Rmd"), 
                  md_document(variant = "markdown_github"))


# 04_promoter_features_profile_plots
rmarkdown::render(file.path("04_promoter_features_profile_plots",
                            "consensus_peak_TSS_meta_plots.Rmd"), 
                  md_document(variant = "markdown_github"))



# 05_promoter_features_lncRNA-vs-mRNA
rmarkdown::render(file.path("05_promoter_features_lncRNA-vs-mRNA",
                            "05_promoter_features_lncRNA_mRNA.Rmd"), 
                  md_document(variant = "markdown_github"))

