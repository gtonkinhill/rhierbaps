#' load_fasta
#'
#' Loads a fasta file into matrix format ready for
#' running the hierBAPS algorithm.
#'
#' @param msa Either the location of a fasta file or ape DNAbin object containing the multiple sequence alignment data to be clustered
#' @param keep.singletons A logical indicating whether to consider singleton mutations in calculating the clusters
#'
#' @return A character matrix with filtered SNP data
#'
#' @examples
#' msa <- system.file("extdata", "seqs.fa", package = "rhierbaps")
#' snp.matrix <- load_fasta(msa)
#' @export
load_fasta <- function(msa, keep.singletons=FALSE) {

  #Check inputs
  if(class(msa)=="character"){
    if (!file.exists(msa)) stop("Invalid msa or the file does not exist!")
    seqs <- ape::read.FASTA(msa)
  } else if(class(msa)=="matrix"){
    seqs <- ape::as.DNAbin(msa)
  } else if(class(msa)=="DNAbin"){
    seqs <- msa
  } else{
    stop("incorrect input for msa!")
  }
  if (!is.logical(keep.singletons)) stop("Invalid keep.singletons! Must be on of TRUE/FALSE.")

  #Load sequences using ape. This does a lot of the checking for us.
  seq_names <- labels(seqs)
  seqs <- as.character(as.matrix(seqs))
  rownames(seqs) <- seq_names
  seqs[is.na(seqs)] <- "-"

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
