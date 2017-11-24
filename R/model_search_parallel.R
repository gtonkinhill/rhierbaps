#' model_search_parallel
#'
#' Clusters DNA alignment using independent loci model
#'
#'
#' @param snp.object A snp.object containing the processed SNP data.
#' @param partition An integer vector indicating an initial starting partition.
#' @param round.types A vector indicating which series of moves to make.
#' @param quiet Whether to suppress progress information (default=FALSE).
#' @param n.extra.rounds The number of additional rounds to perform after the default hierBAPS
#' settings (default=0). If set to Inf it will run until a local optimum is reached
#' (this might take a long time).
#' @param n.cores The number of cores to use.
#'
#' @return
#'
#' @examples
#' snp.matrix <- load_fasta(system.file("extdata", "seqs.fa", package = "rhierbaps"))
#' snp.object <- preproc_alignment(snp.matrix)
#' tmp.hclust <- hclust(as.dist(snp.object$dist), method = 'complete')
#' partition <- cutree(tmp.hclust, k = 20)
#' n.pops <- 20
#' round.types <- c(2*rep(1, n.pops),1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 3, 4,3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3,4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3,4, 3, 4, 3, 4, 3, 4, 3, 4, 1, 1, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4,1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1,2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2,3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3,4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4,1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1,2, 3, 4, 1, 2, 3, 4)
#' model_search_parallel(snp.object, partition, round.types)
#'
model_search_parallel <- function(snp.object, partition, round.types,
                                  quiet=FALSE, n.extra.rounds=0, n.cores=1){
  if(!all(round.types %in% c(1,2,3,4))) stop("Invalid round type!")

  was.updated <- rep(TRUE, 4)
  move.count <- 0
  max.ml <- rhierbaps::calc_log_ml(snp.object, partition)
  if(!quiet){
    cat('\r', paste(c(
      "Round: ", move.count, "/", length(round.types), " Type: ", "none", " Log marginal likelihood: ", max.ml
    ), collapse = ""))
  }

  while(move.count < (length(round.types)+n.extra.rounds)){

    move.count <- move.count + 1
    if(move.count>length(round.types)){
      r <- sample(1:4, 1, replace = TRUE)
    } else {
      r <- round.types[[move.count]]
    }


    n.clusters <- length(unique(partition))
    if (length(unique(partition))<=1){
      next
    }
    if((r==1) && was.updated[[1]]){
      update <- move_units_1(snp.object, partition,
                             n.cores=n.cores)
      if(!update$is.improved){
        was.updated[[1]] <- FALSE
      } else{
        partition <- update$partition
        max.ml <- update$lml
        was.updated <- rep(TRUE, 4)
      }
    } else if(r==2 && was.updated[[2]]){
      update <- join_units_2(snp.object, partition,
                             n.cores=n.cores)
      if(!update$is.improved){
        was.updated[[2]] <- FALSE
      } else{
        partition <- update$partition
        max.ml <- update$lml
        was.updated <- rep(TRUE, 4)
      }
    } else if(r==3 && was.updated[[3]]){
      update <- split_clusters_3(snp.object, partition,
                                 n.cores=n.cores)
      if(!update$is.improved){
        was.updated[[3]] <- FALSE
      } else{
        partition <- update$partition
        max.ml <- update$lml
        was.updated <- rep(TRUE, 4)
      }
    } else if(r==4 && was.updated[[4]]){
      update <- reallocate_units_4(snp.object, partition,
                                   n.cores=n.cores)
      if(!update$is.improved){
        was.updated[[4]] <- FALSE
      } else{
        partition <- update$partition
        max.ml <- update$lml
        was.updated <- rep(TRUE, 4)
      }
    }
    #Print current status
    if(!quiet){
      cat('\r', paste(c(
        "Round: ", move.count, "/", length(round.types), " Type: ", r, " Log marginal likelihood: ", max.ml
      ), collapse = ""))
      flush.console()
    }

    #Check for local convergence
    if (sum(was.updated)==0){
      if (!quiet) print("Converged locally!")
      break
    }
  }
  return(list(partition=partition, lml=rhierbaps::calc_log_ml(snp.object, partition)))
}
