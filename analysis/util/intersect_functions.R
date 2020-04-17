



#Functions we need

# X features_te_total
# Xfeatures overlapping promoters
# >> seperate R file compare gene expression value
#  X features_promoters
# X peak_feature_matrix


#' import peak .bed files as a list
#' 
#' @description 
#' this function will take consensus peak files and name them by the DNA binding protein and return a list
#' 
#' @param consensus_file_path the path to consensus peak files


import_peaks <- function(consensus_file_path = "/Shares/rinn_class/data/CLASS_2020/analysis/01_consensus_peaks") {
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
  
  consensus_peaks <- list()
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
    
    canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")
    for(i in 1:length(peak_list)) {
      peak_list[[i]] <-peak_list[[i]][which(seqnames(peak_list[[i]]) %in% canonical_chr)]
    }
    
    final_peakset <- intersect_peaks(peak_list = peak_list)
    if(length(final_peakset) > 0) {
      final_peakset$name <- paste0(tf, "_", 1:length(final_peakset))
    }
    
    consensus_peaks <- c(consensus_peaks, list(final_peakset))
    names(consensus_peaks)[length(consensus_peaks)] <- tf
  }
  return(consensus_peaks)
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
  genes <- genes[genes$gene_type %in% biotype]
  
  proms <- GenomicRanges::promoters(genes, upstream = upstream, downstream = downstream)
  
  return(proms)
  
}


#' function to subset features for overlapping promomters. 
#' 
#' @description 
#' Take a gencode gtf to subset the biotype of overlapping promoters we want as a set of GRanges
#' 
#' @param gencode_gr
#'  set of genomic features as a GRanges object

get_overlapping_promoters <- function(gencode_gr, upstream = 200, downstream = 0) {
  
  genes <- gencode_gr[gencode_gr$type == "gene"]
  proms <- GenomicRanges::promoters(genes, upstream = upstream, downstream = downstream)
  
  reduced_proms <- GenomicRanges::reduce(proms)
  ov_proms <- GenomicRanges::findOverlaps(proms, reduced_proms, ignore.strand=TRUE)
  ov_proms <- data.frame("gene_promoters_from" = ov_proms@from,
                         "reduced_promoters_to" = ov_proms@to)
  
  ov_proms_summary <- ov_proms %>% 
    group_by(reduced_promoters_to) %>%
    summarize(count = n())
  ov_proms <- merge(ov_proms, ov_proms_summary)
  overlapped <- ov_proms %>% filter(count > 1)
  
  reduced_promoters_overlapping <- reduced_proms[unique(overlapped$reduced_promoters_to)]
  
  overlapped$gene_id <- proms$gene_id[overlapped$gene_promoters_from]
  overlapped$gene_name <- proms$gene_name[overlapped$gene_promoters_from]
  overlapped$gene_type <- proms$gene_type[overlapped$gene_promoters_from]
  overlapped$strand <- strand(proms[overlapped$gene_promoters_from])
  
  overlapped_metadata <- overlapped %>% 
    group_by(reduced_promoters_to, count) %>%
    summarize(gene_id = paste(gene_id, collapse = ";"),
              gene_name = paste(gene_name, collapse = ";"),
              gene_type = paste(gene_type, collapse = ";"),
              strand = paste(as.character(strand), collapse = ";"))
  
  reduced_promoters_overlapping$num_overlaps <- overlapped_metadata$count
  reduced_promoters_overlapping$gene_id <- overlapped_metadata$gene_id
  reduced_promoters_overlapping$gene_name <- overlapped_metadata$gene_name
  reduced_promoters_overlapping$gene_type <- overlapped_metadata$gene_type
  reduced_promoters_overlapping$gene_type <- overlapped_metadata$gene_type
  
  return(reduced_promoters_overlapping) 
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
#' 
#' @param type
#' Return either a matrix of counts over features or a binary occurence matrix

count_peaks_per_feature <- function(features, peak_list, type = "counts") {
  
  if(!(type %in% c("counts", "occurence"))) {
    stop("Type must be either occurence or counts.")
  }
  
  peak_count <- matrix(numeric(), ncol = length(features), nrow = 0)
  
  for(j in 1:length(peak_list)) {
    ov <- countOverlaps(features, peak_list[[j]])
    peak_count <- rbind(peak_count, ov)
    rownames(peak_count)[nrow(peak_count)] <- names(peak_list)[j]
    colnames(peak_count) <- features$gene_id
  }
  
  peak_matrix <- peak_count
  
  if(type == "occurence") {
    peak_occurence <- matrix(as.numeric(peak_count > 0), 
                             nrow = dim(peak_count)[1],
                             ncol = dim(peak_count)[2])
    rownames(peak_occurence) <- rownames(peak_count)
    colnames(peak_occurence) <- colnames(peak_count)
    peak_matrix <- peak_occurence
  }
  
  return(peak_matrix)
  
}





