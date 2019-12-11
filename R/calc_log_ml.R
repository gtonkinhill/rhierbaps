#' calc_log_ml
#'
#' Calculate the log marginal likelihood assuming a Multinomial-Dirichlet distribution
#'
#'
#' @param snp.object A snp.object containing the processed SNP data.
#' @param partition An integer vector indicating a partition of the isolates.
#'
#' @return The log marginal likelihood of the given partition.
#'
calc_log_ml <- function(snp.object, partition){
  #some checks
  if (ncol(snp.object$prior)!=ncol(snp.object$data)) stop("ncol mismatch bwtn prior and data!")
  if (length(partition)!=nrow(snp.object$data)) stop("mismatch bwtn partition and data!")

  #get allele counts for each cluster in partition
  mA <- t(rowsum(1*(snp.object$data=="a"), partition))
  mC <- t(rowsum(1*(snp.object$data=="c"), partition))
  mG <- t(rowsum(1*(snp.object$data=="g"), partition))
  mT <- t(rowsum(1*(snp.object$data=="t"), partition))

  prior <- snp.object$prior
  prior[prior==0] <- 1 #deal with zeros and resulting NAs

  #calculate log marginal likelihood
  term1 <- -lgamma(1 + mA+mC+mG+mT)
  term2 <- lgamma(prior["a", ] + mA) - lgamma(prior["a", ])
  term2 <- term2 + lgamma(prior["c", ] + mC) - lgamma(prior["c", ])
  term2 <- term2 + lgamma(prior["g", ] + mG) - lgamma(prior["g", ])
  term2 <- term2 + lgamma(prior["t", ] + mT) - lgamma(prior["t", ])

  #take sum over all loci and clusters
  ml <- sum(term1 + term2)

  #uniform prior on K:
  #prior prob of each partition is equal so must divide by the number of possible
  #partitions.
  stirling2 <- log_stirling2(nrow(snp.object$data), length(unique(partition)))
  stopifnot(is.finite(stirling2))
  stopifnot(!is.na(stirling2))
  
  ml <- ml - stirling2

  return(ml)
}
