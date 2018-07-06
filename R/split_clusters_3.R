#' split_clusters_3
#'
#' Peform an iteration of the third move in the algorithm. That is split cluster in two
#' and re-allocate one sub-cluster.
#'
#'
#' @param snp.object A snp.object containing the processed SNP data.
#' @param partition An integer vector indicating an initial partition of the isolates.
#' @param threshold The increase in marginal log likelihood required to accept a move.
#' @param min.clust.size Clusters smaller than min.clust.size will not be split.
#' @param n.cores The number of cores to use.
#'
#' @return The best partition after splitting a cluster and re-allocating as well as
#' a boolean value indicating whether a move increased the marginal likelihood.
#'
split_clusters_3 <- function(snp.object, partition, threshold=1e-5,
                             min.clust.size=20, n.cores=1){
  #At the moment this can't create new clusters. This is the same as in the original hierBAPS
  #but it might be worth allowing the creation of new clusters. TODO:Ask Jukka about it.
  return(reallocate_units_4(snp.object, partition, threshold=1e-5,
                            min.clust.size=min.clust.size, split=TRUE,
                            n.cores=n.cores))


}
