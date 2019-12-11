#' preproc_alignment
#'
#' Preprocessed the snp matrix for hierBAPS.
#'
#'
#' @param snp.matrix A matrix containing SNP data. Rows indicate isolates and columns loci.
#'
#' @return an snp.object
#'
preproc_alignment <- function(snp.matrix){
  if(!is.matrix(snp.matrix)) stop("snp.matrix is not a valid matrix")

  n.seq <- nrow(snp.matrix)

  # prior: 1/number of distinct nucleotides observed at that column
  prior <- apply(snp.matrix, 2, function(x) table(factor(x, levels = c("a","c","g","t")),
                                                  exclude = c("-", NA)))

  #TODO: should we be getting rid of unique SNPs here?
  #now ignore conserved columns
  keep <- colSums(prior>0)>1
  
  if(sum(keep)==0){
    #all columns are conserved
    return(NA)
  }
  
  snp.matrix <- snp.matrix[, keep, drop=FALSE]
  prior <- prior[, keep, drop=FALSE]

  #finally generate a matrix for the prior nt values
  prior <- 1*(prior>0)
  prior <- t(t(prior)/colSums(prior))

  #Calculate hamming distance
  orig.dist <- as.matrix(ape::dist.dna(ape::as.DNAbin(snp.matrix),
                                       model = "N", pairwise.deletion = TRUE))

  return(list(
    n.seq = n.seq,
    dist = orig.dist,
    seq.inds = 1:n.seq,
    prior = prior,
    data = snp.matrix
  ))
}
