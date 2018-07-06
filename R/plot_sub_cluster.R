#' plot_sub_cluster
#'
#' Creates a zoom plot using ggtree focusing on a cluster.
#'
#'
#' @param hb.object The resulting object from running hierBAPS
#' @param tree A phylo tree object to plot
#' @param level The level of the subcluster to be considered.
#' @param sub.cluster An integer representing the subcluster to be considered.
#'
#' @examples
#' \donttest{
#' snp.matrix <- load_fasta(system.file("extdata", "seqs.fa", package = "rhierbaps"))
#' newick.file.name <- system.file("extdata", "seqs.fa.treefile", package = "rhierbaps")
#' tree <- phytools::read.newick(newick.file.name)
#' hb.result <- hierBAPS(snp.matrix, max.depth=2, n.pops=20)
#' plot_sub_cluster(hb.result, tree, level = 1, sub.cluster = 9)
#' }
#' @export
plot_sub_cluster <- function(hb.object, tree, level, sub.cluster){
  #Checks
  if ((!is.list(hb.object) || !is.data.frame(hb.object$partition.df)
  ) || !is.list(hb.object$lml.list)) stop("Invalid hb.object!")
  if (!(class(tree)=="phylo")) stop("Invalid tree object!")
  if ((!is.numeric(level)) || (level<1)) stop("Invalid level! Must be a positive integer.")
  if ((!is.numeric(sub.cluster)) || (sub.cluster<1)) stop("Invalid sub.cluster! Must be a positive integer.")


  level <- level+1
  if(!("ggtree" %in%
       rownames(utils::installed.packages()))) stop("This function requires ggtree to be installed")

  cluster.isolate <- hb.object$partition.df$Isolate[hb.object$partition.df[,level]==sub.cluster]

  #Need to create a tempfile to supress the output of gzoom
  ff <- tempfile()
  grDevices::png(filename=ff)
  gg2 <- ggtree::gzoom(tree, focus=which(tree$tip.label %in% cluster.isolate))
  grDevices::dev.off()
  unlink(ff)

  temp_column_id <- paste(c("factor(`level ", level, "`)"), collapse = "")

  p2 <- gg2$p2
  p2 <- ggtree::`%<+%`(p2, hb.object$partition.df)
  p2 <- p2 + ggtree::geom_tippoint(ggplot2::aes_string(color=temp_column_id))
  p2 <- p2+ ggplot2::labs(color=temp_column_id) + ggplot2::theme(legend.position="right") 

  return(ggtree::multiplot(gg2$p1, p2, ncol = 2))
}
