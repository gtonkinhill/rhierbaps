#' preproc_alignment
#'
#' Preprocessed the snp matrix for hierBAPS
#'
#'
#' @param snp.matrix A matrix containing SNP data. Rows indicate isolates and columns loci.
#'
#' @return
#'
#' @examples
#' snp.matrix <- load_fasta(system.file("extdata", "seqs.fa", package = "rhierbaps"))
#' preproc_alignment(snp.matrix)
#'
preproc_alignment <- function(snp.matrix){
  if(class(snp.matrix)!="matrix") stop("snp.matrix is not a valid matrix")

  n.seq <- nrow(snp.matrix)

  # prior: 1/number of distinct nucleotides observed at that column
  prior <- apply(snp.matrix, 2, function(x) table(factor(x, levels = c("a","c","g","t")),
                                                  exclude = c("-", NA)))

  #TODO: should we be getting rid of unique SNPs here?
  #now ignore conserved columns
  keep <- colSums(prior>0)>1
  snp.matrix <- snp.matrix[, keep]
  prior <- prior[, keep]

  #finally generate a matrix for the prior nt values
  prior <- 1*(prior>0)
  prior <- t(t(prior)/colSums(prior))

  #Calculate hamming distance
  orig.dist <- as.matrix(ape::dist.dna(ape::as.DNAbin(snp.matrix),
                                       model = "N"))

  return(list(
    n.seq = n.seq,
    dist = orig.dist,
    seq.inds = 1:n.seq,
    prior = prior,
    data = snp.matrix
  ))
}
