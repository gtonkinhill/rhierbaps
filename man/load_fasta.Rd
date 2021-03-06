% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_fasta.R
\name{load_fasta}
\alias{load_fasta}
\title{load_fasta}
\usage{
load_fasta(msa, keep.singletons = FALSE)
}
\arguments{
\item{msa}{Either the location of a fasta file or ape DNAbin object containing the multiple sequence alignment data to be clustered}

\item{keep.singletons}{A logical indicating whether to consider singleton mutations in calculating the clusters}
}
\value{
A character matrix with filtered SNP data
}
\description{
Loads a fasta file into matrix format ready for
running the hierBAPS algorithm.
}
\examples{
msa <- system.file("extdata", "seqs.fa", package = "rhierbaps")
snp.matrix <- load_fasta(msa)
}
