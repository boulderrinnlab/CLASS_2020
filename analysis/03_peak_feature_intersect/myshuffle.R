my_shuffle <- function (peak_gr, genome) 
{
  chrLens <- genome[names(seqlengths(peak_gr)),"size"]
  nn <- as.vector(seqnames(peak_gr))
  ii <- order(nn)
  w <- width(peak_gr)
  nnt <- table(nn)
  jj <- order(names(nnt))
  nnt <- nnt[jj]
  chrLens <- chrLens[jj]
  ss <- unlist(sapply(1:length(nnt), function(i) sample(chrLens[i], 
                                                        nnt[i])))
  res <- GRanges(seqnames = nn[ii], ranges = IRanges(ss, width = w[ii]), 
                 strand = "*")
  return(res)
}