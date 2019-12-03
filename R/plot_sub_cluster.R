#' plot_sub_cluster
#'
#' Creates a zoom plot using ggtree focusing on a cluster.
#' 
#' @import patchwork
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
  
  # Plot the full tree with the clade highlighted
  hb.object$partition.df$is_in_cluster <- hb.object$partition.df$Isolate %in% cluster.isolate
  full_tree <- ggtree::`%<+%`(ggtree::ggtree(tree), hb.object$partition.df) +
    ggtree::geom_tippoint(ggplot2::aes_string(color="is_in_cluster"), size=0.5) +
    ggplot2::theme(legend.position = "none") + 
    ggplot2::scale_colour_manual(values=c("#000000", "#e31a1c"))
  
  # Subset the tree
  temp_column_id <- paste(c("factor(`level ", level, "`)"), collapse = "")
  sub_tree <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% cluster.isolate])
  sub_tree <- ggtree::`%<+%`(ggtree::ggtree(sub_tree), hb.object$partition.df) +
    ggtree::geom_tippoint(ggplot2::aes_string(color=temp_column_id)) +
    ggplot2::labs(color=temp_column_id) + 
    ggplot2::theme(legend.position="right") +
    ggplot2::scale_color_discrete(name=paste("level ", level, collapse = ""))

  return(full_tree+sub_tree+patchwork::plot_layout(nrow = 1))
}
