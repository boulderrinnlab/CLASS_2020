# Purpose: Trying to make a file with the overlapping binding between the TFs 

library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(ggplot2)
library(tidyverse)


base_path <- "/Shares/rinn_class/data/k562_chip/results/bwa/mergedLibrary/macs/broadPeak/" 

# Making files with peaks for the TFs 

file1 <- rtracklayer::import("/Shares/rinn_class/data/k562_chip/results/bwa/mergedLibrary/macs/broadPeak/ELF1_R1_peaks.broadPeak")
file2 <- rtracklayer::import("/Shares/rinn_class/data/k562_chip/results/bwa/mergedLibrary/macs/broadPeak/ELF1_R2_peaks.broadPeak")

# Find the overlap between the peaks

overlap <- findOverlaps(file1, file2)
overlap


