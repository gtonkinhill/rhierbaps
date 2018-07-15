#' hierBAPS
#'
#' Runs the hierBAPS algorithm of Cheng et al. 2013
#'
#'
#' @param snp.matrix Character matrix of aligned sequences produced by \link{load_fasta}.
#' @param max.depth Maximum depth of hierarchical search (default = 2).
#' @param n.pops Maximum number of populations in the data (default = number of isolates/5)
#' @param quiet Whether to suppress progress information (default=FALSE).
#' @param n.extra.rounds The number of additional rounds to perform after the default hierBAPS
#' settings (default=0). If set to Inf it will run until a local optimum is reached
#' (this might take a long time).
#' @param n.cores The number of cores to use.
#'
#' @return A list containing a dataframe indicating an assignment of each sequence
#' to hierarchical clusters as well as the log marginal likelihoods for each level.
#'
#' @examples
#' snp.matrix <- load_fasta(system.file("extdata", "small_seqs.fa", package = "rhierbaps"))
#' hb <- hierBAPS(snp.matrix, max.depth=2, n.pops=20, quiet=FALSE)
#'
#' \donttest{
#' snp.matrix <- load_fasta(system.file("extdata", "seqs.fa", package = "rhierbaps"))
#' system.time({hb <- hierBAPS(snp.matrix, max.depth=2, n.pops=20, quiet=FALSE)})
#' }
#'
#'@author Gerry Tonkin-Hill
#'@references Cheng, Lu, Thomas R. Connor, Jukka Sirén, David M. Aanensen, and Jukka Corander. 2013. “Hierarchical and Spatially Explicit Clustering of DNA Sequences with BAPS Software.” Molecular Biology and Evolution 30 (5): 1224–28.
#'
#' @export
hierBAPS <- function(snp.matrix, max.depth=2, n.pops=floor(nrow(snp.matrix)/5),
                     quiet=FALSE, n.extra.rounds=0, n.cores=1){

  #Check inputs
  if (class(snp.matrix)!="matrix") stop("snp.matrix is not a matrix!")
  if ((!is.numeric(max.depth)) || (max.depth<1)) stop("Invalid max.depth! Must be a positive integer.")
  if ((!is.numeric(n.pops)) || (n.pops<1)) stop("Invalid n.pops! Must be a positive integer.")
  if (!is.logical(quiet)) stop("Invalid quiet! Must be one of TRUE/FALSE.")
  if ((!is.numeric(n.extra.rounds)) || (n.extra.rounds<0)) stop("Invalid n.extra.rounds!
                                                                Must be a non-negative integer.")
  if ((!is.numeric(n.cores)) || (n.cores<1)) stop("Invalid n.cores! Must be a positive integer.")

  # search operators
  round.types <- c(2*rep(1, n.pops),
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 3, 4,
                 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3,
                 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4,
                 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3,
                 4, 3, 4, 3, 4, 3, 4, 3, 4, 1, 1, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4,
                 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1,
                 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2,
                 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3,
                 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4,
                 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1,
                 2, 3, 4, 1, 2, 3, 4)
  
  #Don't split clusters with the less than 4 memebers
  MIN.CLUSTER.SIZE <- 4

  #iterate over levels
  all.partition.matrix <- matrix(0, nrow=nrow(snp.matrix), ncol=max.depth)
  lml.list <- list()

  for (cur.depth in seq(0, max.depth-1)){

    if (cur.depth==0){
      snp.object <- preproc_alignment(snp.matrix)
      if(is.na(snp.object[[1]])){
        stop("All sites are conserved!")
      }
      cur.part <- rep(1, snp.object$n.seq)
      avail.cluster.ids <- 1
      snp.object$heds <- rownames(snp.matrix)
    } else{
      cur.part <- all.partition.matrix[snp.object$seq.inds, cur.depth]
      avail.cluster.ids <- unique(cur.part)
    }

    if(!quiet) print(paste(c("---- Current depth: ", cur.depth, " ----"), collapse = ""))

    if(length(cur.part)==0){
      if(!quiet) print(paste(c("all sequences are clustered. Quit! There are ", cur.depth,
                  " layers in total."), collapse = ""))
      break
    }

    local.label.offset <- 1

    for (i in 1:length(avail.cluster.ids)){
      cluid <- avail.cluster.ids[i]
      if(!quiet) print(paste(c("Current depth: ", cur.depth, " Cluster ID: ", cluid), collapse = ""))

      if (sum(cur.part==cluid) < MIN.CLUSTER.SIZE){
        #we dont split any further.
        all.partition.matrix[cur.part==cluid, cur.depth+1] <- local.label.offset
        local.label.offset <- local.label.offset+1
        if(length(lml.list)<(cur.depth+1)){
          lml.list[[cur.depth+1]] <- NA
        } else {
          lml.list[[cur.depth+1]] <- c(lml.list[[cur.depth+1]], NA)
        }
        next
      }

      if (cur.depth==0){
        tmp.snp.object <- snp.object
      } else {
        tmp.snp.object <- preproc_alignment(snp.object$data[cur.part==cluid, ])
        if(is.na(tmp.snp.object[[1]])){
          #all sites are conserved at this partition so dont split further
          all.partition.matrix[cur.part==cluid, cur.depth+1] <- local.label.offset
          local.label.offset <- local.label.offset+1
          if(length(lml.list)<(cur.depth+1)){
            lml.list[[cur.depth+1]] <- NA
          } else {
            lml.list[[cur.depth+1]] <- c(lml.list[[cur.depth+1]], NA)
          }
          next
        }
      }
      
      

      tmp.z.hclust <- stats::hclust(stats::as.dist(tmp.snp.object$dist), method = 'complete')

      if (tmp.snp.object$n.seq > (3*n.pops)){
        tmp.init.part <- stats::cutree(tmp.z.hclust, k = n.pops)
      } else {
        tmp.num = min(floor(tmp.snp.object$n.seq/2), n.pops)
        tmp.init.part <- stats::cutree(tmp.z.hclust, k = tmp.num)
      }

      temp_partition = model_search_parallel(tmp.snp.object, tmp.init.part, round.types,
                                             quiet, n.extra.rounds, n.cores)

      if (!quiet){
        print(paste(c("Best partition: Nclusters ", length(unique(temp_partition$partition)),
                      " Log(ml*prior) ", temp_partition$lml), collapse = ""))
      }

      if(length(lml.list)<(cur.depth+1)){
        lml.list[[cur.depth+1]] <- temp_partition$lml
      } else {
        lml.list[[cur.depth+1]] <- c(lml.list[[cur.depth+1]], temp_partition$lml)
      }


      for (clust in unique(temp_partition$partition)){
        cluster_index <- rep(FALSE, length(cur.part))
        cluster_index[cur.part==cluid][clust==temp_partition$partition] <- TRUE
        all.partition.matrix[cluster_index, cur.depth+1] <- local.label.offset
        local.label.offset <- local.label.offset+1
      }
    }

    names(lml.list[[cur.depth+1]]) <- avail.cluster.ids

  }
  names(lml.list) <- paste("Depth", seq(0, max.depth-1))
  partition.df <- data.frame(isolates=snp.object$heds, all.partition.matrix)
  colnames(partition.df) <- c("Isolate", paste("level", 1:ncol(all.partition.matrix)))
  return(list(partition.df=partition.df,
              lml.list=lml.list))
}
