#' log_stirling2
#'
#' Calculate the log of the Stirling number of the second kind
#' Adapted from the function 'Strlng2' in the CryptRndTest package
#'
#' @param n positive integer greater than zero
#' @param k positive integer between 1 and n
#'
#' @return The log of the Stirling number of the second kind.
#' 
#' @references Bleick, W.W., Wang, P.C.C., Asymptotics of Stirling Numbers of the Second Kind. Proceedings of the American Mathematical Society (1974), 42(2), 575–580.
#'
log_stirling2 <- function(n, k){
  #some checks
  if (!is.numeric(n) || n<0) stop("n must be a positive integer!")
  if (!is.numeric(k) || (k<0 || k>n)) stop("k must be a positive integer between 1 and n!")
  
  if (k == 1) {
    topl = 1
  }
  else {
    z = 1:n
    nu = n/k
    r = 1:(n - k)
    G = -LambertW::W(-nu * exp(-nu))
    topl = log(sqrt(n - k)) + (n - k) * log((n - k)/exp(1)) + 
      (sum(log(z)) - sum(log(r)) - sum(log(1:k))) - log(sqrt(n * 
                                                               (1 - G))) - k * log(G) - (n - k) * log(nu - G)
  }
  
  return(topl)
}
