#' log_stirling2
#'
#' @param n number of objects
#' @param k number of partitions
#'
#' @return log of the Stirling number of the second kind
#' 
#'
log_stirling2 <- function(n, k){
  if(!is.numeric(n)) stop("n is not numeric!")
  if(!is.numeric(k)) stop("k is not numeric!")
  if(k>n) stop("k must be less than n!")
  
  v <- n/k
  G <- lambertW(-v*exp(-v))
  
  lS2 <- log(sqrt((v-1)/(v*(1-G)))) +
    (n-k)*(log(v-1)-log(v-G)) +
    n*log(k)-k*log(n) +
    k*(1-G) +
    lchoose(n, k)

  return(lS2)
}

# This function was written by Ben Bolker and taken from https://stat.ethz.ch/pipermail/r-help/2003-November/042793.html
lambertW = function(z,b=0,maxiter=10,eps=.Machine$double.eps,
                    min.imag=1e-9) {
  if (any(round(Re(b)) != b))
    stop("branch number for W must be an integer")
  if (!is.complex(z) && any(z<0)) z=as.complex(z)
  ## series expansion about -1/e
  ##
  ## p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  ## w = (11/72)*p;
  ## w = (w - 1/3).*p;
  ## w = (w + 1).*p - 1
  ##
  ## first-order version suffices:
  ##
  w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
  ## asymptotic expansion at 0 and Inf
  ##
  v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
  v = v - log(v + as.numeric(v==0))
  ## choose strategy for initial guess
  ##
  c = abs(z + exp(-1));
  c = (c > 1.45 - 1.1*abs(b));
  c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
  w = (1 - c)*w + c*v
  ## Halley iteration
  ##
  for (n in 1:maxiter) {
    p = exp(w)
    t = w*p - z
    f = (w != -1)
    t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
    w = w - t
    if (abs(Re(t)) < (2.48*eps)*(1.0 + abs(Re(w)))
        && abs(Im(t)) < (2.48*eps)*(1.0 + abs(Im(w))))
      break
  }
  if (n==maxiter) warning(paste("iteration limit (",maxiter,
                                ") reached, result of W may be inaccurate",sep=""))
  if (all(Im(w)<min.imag)) w = as.numeric(w)
  return(w)
}
