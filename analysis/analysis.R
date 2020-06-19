######
# This script will run the full analysis using the singularity container 
#####

#### SETUP
library(rmarkdown)
# TODO: potentially make all of these into readme.md files
# so that they render automatically in github when the file path is
# clicked. Alternatively, render all the markdowns into a github.io site
# according to these instructions:
# https://nicolas-van.github.io/easy-markdown-to-github-pages/

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


# 06_repeat_features
rmarkdown::render(file.path("06_repeat_features",
                            "06_repeat_features.Rmd"), 
                  md_document(variant = "markdown_github"))


# 07_binding_vs_expression
rmarkdown::render(file.path("07_binding_versus_expression",
                            "07_binding_vs_expression.Rmd"), 
                  md_document(variant = "markdown_github"))


# 08_defining_reservoirs
rmarkdown::render(file.path("08_defining_reservoirs",
                            "08_defining_reservoirs.Rmd"), 
                  md_document(variant = "markdown_github"))


# 09_reservoir_binding_properties
rmarkdown::render(file.path("09_reservoir_binding_properties",
                            "09_reservoir_binding_properties.Rmd"), 
                  md_document(variant = "markdown_github"))


# 10_reservoir_nascent_txn
rmarkdown::render(file.path("10_reservoir_nascent_txn",
                            "10_reservoir_nascent_txt.Rmd"), 
                  md_document(variant = "markdown_github"))


# 03_umap_with_metadata
rmarkdown::render(file.path("03_global_clustering",
                            "umap_with_metadata.Rmd"), 
                  md_document(variant = "markdown_github"))

