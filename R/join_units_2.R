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
#' @param comb.chache a matrix recording previous marginal llks of combining clusters
#'
#' @return The best partition after combining two clusters as well as
#' a boolean value indicating whether a move increased the marginal likelihood.
#'
join_units_2 <- function(snp.object, partition, threshold=1e-5, n.cores=1, comb.chache=NULL){

  #some checks
  if (ncol(snp.object$prior)!=ncol(snp.object$data)) stop("ncol mismatch bwtn prior and data!")
  if (length(partition)!=nrow(snp.object$data)) stop("mismatch bwtn partition and data!")

  max_ml <- calc_log_ml(snp.object, partition)
  is.improved <- FALSE
  clusters <- unique(partition)

  stopifnot(length(clusters) > 1)
  
  
  combinations <- utils::combn(clusters, 2, simplify = TRUE)
  
  if(is.null(comb.chache)){
    comb.chache <- matrix(NA, nrow=max(clusters), ncol=max(clusters))
  }
  
  temp.combinations <- combinations[, is.na(comb.chache[t(combinations)]), drop=FALSE]
  
  temp_mls <- parallel::mcmapply(function(p1, p2){
    temp_partition <- partition
    temp_partition[temp_partition==p2] <- p1
    return(calc_log_ml(snp.object, temp_partition))
  }, temp.combinations[1,], temp.combinations[2,], mc.cores = n.cores)
  
  comb.chache[t(temp.combinations)] <- temp_mls

  arg.max <- which(comb.chache == max(comb.chache, na.rm = TRUE), arr.ind = TRUE)[1,]

  if(comb.chache[arg.max[[1]], arg.max[[2]]] > (max_ml+threshold)){
    partition[partition==arg.max[[2]]] <- arg.max[[1]]
    is.improved <- TRUE
    if (length(clusters)<=2){
      diff <- (comb.chache[arg.max[[1]], arg.max[[2]]] - 
                 max_ml + 
                 2*log_stirling2(nrow(snp.object$data), length(clusters)-1) - 
                 log_stirling2(nrow(snp.object$data), length(clusters))
      )
    } else {
      diff <- (comb.chache[arg.max[[1]], arg.max[[2]]] - 
                 max_ml + 
                 2*log_stirling2(nrow(snp.object$data), length(clusters)-1) - 
                 log_stirling2(nrow(snp.object$data), length(clusters)) -
                 log_stirling2(nrow(snp.object$data), length(clusters)-2)
      )
    }
    
    max_ml <- comb.chache[arg.max[[1]], arg.max[[2]]]
    comb.chache <- comb.chache + diff
    comb.chache[arg.max[[1]], ] <- NA
    comb.chache[arg.max[[2]], ] <- NA
    comb.chache[, arg.max[[1]]] <- NA
    comb.chache[, arg.max[[2]]] <- NA
  }

  return(list(partition=partition, is.improved=is.improved, lml=max_ml, comb.chache=comb.chache))
}
