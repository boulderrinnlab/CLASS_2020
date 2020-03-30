# We will use this to practice debugging.

count_peaks_per_feature_r <- function(features, peak_list) {
  
  peak_count <- matrix(numeric(), ncol = length(features), nrow = 0)
  
  for(j in 1:length(peak_list)) {
    ov <- countOverlaps(features, peak_list[[j]])
    peak_count <- rbind(peak_count, ov)
    rownames(peak_count)[nrow(peak_count)] <- names(peak_list)[j]
    colnames(peak_count) <- features$gene_id
  }
  return(peak_count)
}
