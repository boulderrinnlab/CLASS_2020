



#Functions we need
# features_genebody_3kb
# features_te_total
#  X features_promoters
# X peak_feature_matrix


#' import peak .bed files as a list
#' 
#' @description 
#' this function will take consensus peak files and name them by the DNA binding protein and return a list
#' 
#' @param consensus_file_path the path to consensus peak files
#' 

import_peaks <- function(consensus_file_path = "/Shares/rinn_class/data/k562_chip/analysis/00_consensus_peaks/results/") {
  peak_files <- list.files(consensus_file_path, full.names = T)
  file_names <- str_extract(peak_files, "[\\w-]+\\.bed")
  tf_name <- str_extract(file_names, "^[^_]+(?=_)")
  
  
  peak_list <- c()
  for(i in 1:length(peak_files)) {
    # Import peaks
    peaks <- rtracklayer::import(peak_files[i])
    # Append this GRanges object to the of the list.
    peak_list <- c(peak_list, peaks)
    # Name the list elements by their TF name.
    names(peak_list)[length(peak_list)] <- tf_name[i]
  }
  return(peak_list)
}



#' intersect replicates into a "consensus peak list" 
#' 
#' @description 
#' this function will take the intersect and union of peak widths across replicates for a given DNA binding protein. the function that will take a list of granges objects and return 
# one granges object with merged peaks that are in all replicates
#' 
#' @param 
#'  the path to consensus peak files
#' # We're going to iterate over all the files to make it work. 

create_consensus_peaks <- function(broadpeakfilepath = "/Shares/rinn_class/data/k562_chip/results/bwa/mergedLibrary/macs/broadPeak") {
  
  
  fl <- list.files(broadpeakfilepath, 
                   full.names=TRUE)
  fl <- fl[grep("peaks.broadPeak", fl)]
  
  tf_name <- sapply(fl, function(x){
    y <-  unlist(strsplit(x, "/"))[[11]]
    unlist(strsplit(y, "_"))[[1]]
  })
  
  unique_tf <- unique(tf_name)
  
  # This for loop will iterate over all dna binding proteins.
  for(i in 1:length(unique_tf)) {
    
    # load all the peak files corresponding to this dna binding proteins.
    tf <- unique_tf[i]
    tf_index <- grep(tf, tf_name)
    tf_files <- fl[tf_index]
    
    peak_list <- c()
    for(j in 1:length(tf_files)) {
      # See the read peaks function to know what subfunctions are called.
      peak_list <- c(peak_list, read_peaks(tf_files[j]))
    }
    
    final_peakset <- intersect_peaks(peak_list = peak_list)
    if(length(final_peakset) > 0) {
      final_peakset$name <- paste0(tf, "_", 1:length(final_peakset))
    }
    # write out that peakset as a bed file. 
    #rtracklayer::export(final_peakset, paste0("results/", tf, "_consensus_peaks.bed"))
  }
  return(final_peakset)
}

# TODO: refactor
read_peaks <- function(broad_peak_file) {
  # A broad peak file is just a tab separated file 
  dat <- read.table(broad_peak_file, sep = "\t")
  gr <- GRanges(seqnames = dat$V1,
                ranges = IRanges(start=dat$V2,end=dat$V3))
  return(gr)
}
# This is the function that will be doing the core of the work here. 
# When two peaks intercept, we will take their outer boundaries to be the new
# peak -- using the reduce function.
intersect_peaks <- function(peak_list) {
  combined_peaks <- peak_list[[1]]
  for(i in 2:length(peak_list)) {
    pl_ov <- findOverlaps(combined_peaks, peak_list[[i]])
    pl1 <- combined_peaks[unique(pl_ov@from)]
    pl2 <- peak_list[[i]][unique(pl_ov@to)]
    combined_peaks <- GenomicRanges::reduce(union(pl1, pl2))
  }
  return(combined_peaks)
}



#' function finds overlaps between consensus_peaks and genomic features 
#' 
#' @description get_overlapping_peaks
#' this function will intersect the consesus_peaks with gene_body, promoters, mRNA_promoters, lncRNA_promoters, te_family
#' 
#' @param features
#'  set of genomic features as a GRanges object
#'  
#' @param peak_list
#' #list of peaks of dna binding proteins that will be intersected

get_overlapping_peaks <- function(features, peak_list){
  
  feature_peaks <- c()
  
  overlaps_list <- c()
  for(j in 1:length(peak_list)) {
    ov <- findOverlaps(peak_list[[j]], features)
    overlapping_peaks <- peak_list[[j]][unique(ov@from)]
    overlaps_list <- c(overlaps_list, overlapping_peaks)
    names(overlaps_list)[length(overlaps_list)] <- names(peak_list)[j]
  }
  
  feature_peaks <- c(feature_peaks, list(overlaps_list))
  names(feature_peaks)[length(feature_peaks)] <- names(feature_sets)[i]
  
  return(feature_peaks)
} 


# subset rmsk

import_repeatmasker <- function(rmsk_file = "/Shares/rinn_class/data/k562_chip/util/rmsk.txt") {
  
  rmsk <- read.table(file = rmsk_file)
  colnames(rmsk) <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                      "genoName", "genoStart", "genoEnd", "genoLeft", "strand",
                      "repName", "repClass", "repFamily", "repStart",	"repEnd",
                      "repLeft",	"id")
  
  rmsk_gr <- GRanges(seqnames = rmsk$genoName,
                     IRanges(rmsk$genoStart, rmsk$genoEnd),
                     strand = rmsk$strand)
  
  # Add metadata colums
  rmsk_gr$rep_class <- rmsk$repClass
  rmsk_gr$rep_family <- rmsk$repFamily
  rmsk_gr$rep_name <- rmsk$repName
  
  return(rmsk_gr)
}

subset_rmsk <- function(rmsk_gr, rep_level = "family") {
  # rep_level needs to be one of "class", "family", or "name"
  if(!(rep_level %in% c("class", "family", "name"))) {
    stop("Repeat level needs to be either: class, family, or name.")
  } 
  
  level_column <- grep(rep_level, names(rmsk_gr@elementMetadata))
  rep_levels <- unique(rmsk_gr@elementMetadata[,level_column])
  
  rmsk_list <- c()
  for(i in 1:length(rep_levels)) {
    rmsk_list <- c(rmsk_list, list(rmsk_gr[rmsk_gr@elementMetadata[,level_column] == rep_levels[i]]))
    names(rmsk_list)[length(rmsk_list)] <- rep_levels[[i]]
  }
  return(rmsk_list)
}













#' function to subset features for promomters. 
#' 
#' @description feature_subset
#' Take a gencode gtf to subset the biotype of promoters we want as a set of GRanges
#' 
#' @param gencode_gr
#'  set of genomic features as a GRanges object
#'  
#' @param biotype
#' this takes "lncRNA" or "protein-coding" as input for promoter type
#'
#' @param upstream
#'To add upstream sequence to feature subset
#'
#' @param downstream
#'To add downstream sequence to feature subset

get_promoter_regions <- function(gencode_gr, biotype, upstream = 3e3, downstream = 3e3) {
  
  genes <- gencode_gr[gencode_gr$type == "gene"]
  genes <- genes[which(genes$gene_type == biotype)]
  
  proms <- GenomicRanges::promoters(genes, upstream = upstream, downstream = downstream)
  
  return(proms)
  
}


#' function to summarize the number of events in features on each individual promoter. 
#' 
#' @description 
#' Take a gencode gtf to subset the biotype of promoters we want as a set of GRanges
#' 
#' @param features
#' set of genomic features as a GRanges object
#'  
#' @param peak_list
#' #list of peaks of dna binding proteins that will be intersected


count_peaks_per_feature <- function(features, peak_list) {
  
  peak_count_list <- list()
  
  peak_count <- matrix(numeric(), ncol = length(features), nrow = 0)
  
  for(j in 1:length(peak_list)) {
    ov <- countOverlaps(features, peak_list[[j]])
    peak_count <- rbind(peak_count, ov)
    rownames(peak_count)[nrow(peak_count)] <- names(peak_list)[j]
    colnames(peak_count) <- features$gene_id
  }
  
  peak_count_list <- c(peak_count_list, list(peak_count))
  names(peak_count_list)[i] <- names(feature_sets)[i]
  
  return(peak_count_list)
  
}



