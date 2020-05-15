######
# This script will run the full analysis using the singularity container 
#####

#### SETUP
library(rmarkdown)


# 01_consensus_peaks
rmarkdown::render("01_consensus_peaks/consensus_peak_set.Rmd", md_document(variant = "markdown_github"))

# 02_tss_metaplots
rmarkdown::render("02_tss_metaplots/tss_metaplots.Rmd", md_document(variant = "markdown_github"))

# 03_clustering
rmarkdown::render("03_clustering/global_clustering.Rmd", md_document(variant = "markdown_github"))

# 04_permutation_test
rmarkdown::render("04_permutation_test/permutation_test_by_feature_type.Rmd", md_document(variant = "markdown_github"))

# 05_rnaseq_expression
rmarkdown::render("05_rnaseq_expression/expression_vs_binding.Rmd", md_document(variant = "markdown_github"))

# 06_nascent_expression
rmarkdown::render("06_nascent_expression/nascent_expression.Rmd", md_document(variant = "markdown_github"))

# 07_reservoir_properties
rmarkdown::render("07_reservoir_properties/factors_on_reservoirs.Rmd", md_document(variant = "markdown_github"))
