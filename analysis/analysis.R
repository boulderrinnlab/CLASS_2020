######
# This script will run the full analysis using the singularity container 
#####

#### SETUP
library(rmarkdown)


# 01_consensus_peaks
rmarkdown::render("01_consensus_peaks/01_consensus_peaks.Rmd", md_document(variant = "markdown_github"))
