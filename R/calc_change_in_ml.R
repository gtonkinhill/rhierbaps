#' calc_change_in_ml
#'
#' Calculate the change in the log marginal likelihood after moving index to each possible cluster
#'
#'
#' @param snp.object A snp.object containing the processed SNP data.
#' @param partition An integer vector indicating a partition of the isolates.
#' @param indexes Indexes of the isolates to be moved (must come from one cluster.)
#'
#' @return the best cluster to move indexes to.
#'
calc_change_in_ml <- function(snp.object, partition, indexes){
  #some checks
  if (ncol(snp.object$prior)!=ncol(snp.object$data)) stop("ncol mismatch bwtn prior and data!")
  if (length(partition)!=nrow(snp.object$data)) stop("mismatch bwtn partition and data!")
  if (!(all(indexes %in% 1:length(partition)))) stop("indexes outside of partiton range!")
  if (!(class(partition)=="integer")) stop("partition is not an integer vector!")

  original_cluster <- unique(partition[indexes])
  if(length(original_cluster)!=1) stop("there was not a unique cluster in the index set!")

  #create temporary partition with indexes in a seperate cluster
  temp_partition <- partition
  temp_cluster <- max(partition)+1
  temp_partition[indexes] <- temp_cluster

  #get allele counts for each cluster in partition
  mA <- t(rowsum(1*(snp.object$data=="a"), temp_partition))
  mC <- t(rowsum(1*(snp.object$data=="c"), temp_partition))
  mG <- t(rowsum(1*(snp.object$data=="g"), temp_partition))
  mT <- t(rowsum(1*(snp.object$data=="t"), temp_partition))

  #now add counts from the temp cluster to every other cluster and take away from its own.
  mA <- mA[,colnames(mA)!=temp_cluster, drop=FALSE] + mA[, colnames(mA)==temp_cluster]
  mC <- mC[,colnames(mC)!=temp_cluster, drop=FALSE] + mC[, colnames(mC)==temp_cluster]
  mG <- mG[,colnames(mG)!=temp_cluster, drop=FALSE] + mG[, colnames(mG)==temp_cluster]
  mT <- mT[,colnames(mT)!=temp_cluster, drop=FALSE] + mT[, colnames(mT)==temp_cluster]

  prior <- snp.object$prior
  prior[prior==0] <- 1 #deal with zeros and resulting NAs

  #calculate log marginal likelihood
  term1 <- -lgamma(1 + mA+mC+mG+mT)
  term2 <- lgamma(prior["a", ] + mA) - lgamma(prior["a", ])
  term2 <- term2 + lgamma(prior["c", ] + mC) - lgamma(prior["c", ])
  term2 <- term2 + lgamma(prior["g", ] + mG) - lgamma(prior["g", ])
  term2 <- term2 + lgamma(prior["t", ] + mT) - lgamma(prior["t", ])

  new_columnwise_lml <- colSums(term1 + term2)
  names(new_columnwise_lml) <- colnames(mA)

  #now calculate original columnwise llm before the move
  #get allele counts for each cluster in partition
  mA <- t(rowsum(1*(snp.object$data=="a"), partition))
  mC <- t(rowsum(1*(snp.object$data=="c"), partition))
  mG <- t(rowsum(1*(snp.object$data=="g"), partition))
  mT <- t(rowsum(1*(snp.object$data=="t"), partition))
  term1 <- -lgamma(1 + mA+mC+mG+mT)
  term2 <- lgamma(prior["a", ] + mA) - lgamma(prior["a", ])
  term2 <- term2 + lgamma(prior["c", ] + mC) - lgamma(prior["c", ])
  term2 <- term2 + lgamma(prior["g", ] + mG) - lgamma(prior["g", ])
  term2 <- term2 + lgamma(prior["t", ] + mT) - lgamma(prior["t", ])
  orignal_columnwise_lml <- colSums(term1 + term2)
  names(orignal_columnwise_lml) <- colnames(mA)
  new_columnwise_lml <- new_columnwise_lml[names(orignal_columnwise_lml)]

  diff_columnwise_lml <- orignal_columnwise_lml - new_columnwise_lml
  diff_columnwise_lml[original_cluster] <- Inf

  best_column <- names(orignal_columnwise_lml)[which.min(diff_columnwise_lml)]

  return(best_column)
}
