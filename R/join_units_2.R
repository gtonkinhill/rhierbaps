#' join_units_2
#'
#' Peform an iteration of the second move in the algorithm. That is combine two clusters
#' to improve the marginal likelihood.
#'
#'
#' @param snp.object A snp.object containing the processed SNP data.
#' @param partition An integer vector indicating an initial partition of the isolates.
#' @param threshold The increase in marginal log likelihood required to accept a move.
#' @param n.cores The number of cores to use.
#'
#' @return The best partition after combining two clusters as well as
#' a boolean value indicating whether a move increased the marginal likelihood.
#'
#' @examples
#' snp.matrix <- load_fasta(system.file("extdata", "seqs.fa", package = "rhierbaps"))
#' snp.object <- preproc_alignment(snp.matrix)
#' tmp.hclust <- hclust(as.dist(snp.object$dist), method = 'complete')
#' partition <- cutree(tmp.hclust, k = 20)
#' rhierbaps:::join_units_2(snp.object, partition)
join_units_2 <- function(snp.object, partition, threshold=1e-5, n.cores=1){

  #some checks
  if (ncol(snp.object$prior)!=ncol(snp.object$data)) stop("ncol mismatch bwtn prior and data!")
  if (length(partition)!=nrow(snp.object$data)) stop("mismatch bwtn partition and data!")

  max_ml <- rhierbaps::calc_log_ml(snp.object, partition)
  is.improved <- FALSE
  clusters <- unique(partition)

  stopifnot(length(clusters) > 1)

  combinations <- utils::combn(clusters, 2, simplify = TRUE)

  temp_mls <- parallel::mcmapply(function(p1, p2){
    temp_partition <- partition
    temp_partition[temp_partition==p2] <- p1
    return(calc_log_ml(snp.object, temp_partition))
  }, combinations[1,], combinations[2,], mc.cores = n.cores)

  arg.max <- which.max(temp_mls)

  if(temp_mls[arg.max] > (max_ml+threshold)){
    pair <- combinations[, arg.max]
    partition[partition==pair[[2]]] <- pair[[1]]
    is.improved <- TRUE
    max_ml <- temp_mls[arg.max]
  }

  return(list(partition=partition, is.improved=is.improved, lml=max_ml))
}
