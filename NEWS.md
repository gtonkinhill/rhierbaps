# rhierbaps 1.1.2
Fixed issue with stirling2 dependency causing the algorithm to stop prematurely for very large alignments.

# rhierbaps 1.1.1
Removed call to ggtree in introduction vignette as waiting on patch for ggtree given the new version of tibble.

# rhierbaps 1.1.0
* minor bug fix that caused an error when a cluster contained only a single SNP
* added functionality to load in ape DNAbin objects
* added functionality to calculate individual assignment probabilities

# Initial release 1.0.0devtools::release()