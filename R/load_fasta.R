#' load_fasta
#'
#' Loads a fasta file into matrix format ready for
#' running the hierBAPS algorithm.
#'
#' @param fasta.file A fasta file containing the sequence data to be clustered
#' @param keep.singletons A logical indicating whether to consider singleton mutations in calculating the clusters
#'
#' @return A character matrix with filtered SNP data
#'
#' @examples
#' fasta.file <- system.file("extdata", "seqs.fa", package = "rhierbaps")
#' snp.matrix <- load_fasta(fasta.file)
#' @export
load_fasta <- function(fasta.file, keep.singletons=FALSE) {

  #Check inputs
  if (!file.exists(fasta.file)) stop("Invalid fasta.file or the file does not exist!")
  if (!is.logical(keep.singletons)) stop("Invalid keep.singletons! Must be on of TRUE/FALSE.")

  #Load sequences using ape. This does a lot of the checking for us.
  seqs <- ape::read.FASTA(fasta.file)
  seq_names <- labels(seqs)
  seqs <- as.character(as.matrix(seqs))
  rownames(seqs) <- seq_names

  if (nrow(seqs)<3) stop("Less than 3 sequences!")
  warning("Characters not in acgtnACGTN- will be treated as missing (-)...")

  #Remove conserved columns
  conserved <- colSums(t(t(seqs)==seqs[1,]))==nrow(seqs)
  seqs <- seqs[, !conserved]

  if(!keep.singletons){
    #remove singletons as they are uninformative in the algorithm
    is_singleton <- apply(seqs, 2, function(x){
      tab <- table(x)
      return(x %in% names(tab)[tab==1])
    })
    seqs[is_singleton] <- "-"
  }

  #Convert gaps and unknowns to same symbol
  seqs[seqs=="n"] <- "-"

  return(seqs)
}
