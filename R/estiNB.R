# copyright: Xi Wang (xi.wang@newcastle.edu.au)
# two subfunctions to solve D-equation 
phifun <- function(x, ctable, mu) {
  n_exon <- nrow(ctable)
  n_sample <- ncol(ctable)
  f <- (ctable + x) * log((ctable + x)/(mu + x)) - 
    ctable * log(ctable/mu) 
  sum(f) + (n_exon-1)*(n_sample-1)/2
}

dphifun <- function(x, ctable, mu) {
  f <- log((ctable + x)/(mu + x)) + (mu - ctable) / (mu + x)
  sum(f)
}

estiPhi <- function(ctable) {
  n_exon <- nrow(ctable)
  n_sample <- ncol(ctable)
  m <- colSums(ctable)
	
	# initialize w and prob 
	w <- m/sum(m)
	prob <- rowSums(ctable) / sum(m)
  # initialize dispersion 
	x <- 1e-6    ###initial value 1e-6, x=1/phi
  mu <- prob %*% t(m)  
  xnew <- x - phifun(x, ctable, mu)/dphifun(x, ctable, mu)
	phi <- 1/xnew
	
	if(phi == 0) {
		iternum <- 1
		obj <- list(phi, w, prob, iternum)
		names(obj) <- c("phi", "weights", "prob", "iternum")
    return(obj)
	} 
  
  iter_max <- 1000
  iternum_phi <- 1
  # estimate phi using Newton iteration method
  while (iternum_phi <= iter_max) {
    x <- xnew
    phi_fun <- phifun(x, ctable, mu)
    if( abs(phi_fun) < 1e-8 ) {break} 
    iternum_phi <- iternum_phi + 1
    xnew <- x - phi_fun / dphifun(x, ctable, mu)
    if ( abs(xnew/x - 1) < 1e-6 ) {break} 
    if(xnew > 1e10) { xnew <- Inf; break }
    if(xnew < -1e10) { xnew <- -Inf; break } 
  }
  phinew <- 1/xnew
  
  # iterate to get phi and prob
  iternum <- 1
  while(iternum <= iter_max) {
    iternum <- iternum + 1
    phi <- phinew
    # update w 
    p2 <- sum(prob^2)
    w <- m / (1 + phi * p2 * m )
    w <- w/sum(w)
    # update prob
    prob <- ctable %*% (w / m)
    # update phi
    x <- 1e-6    ###initial value 1e-6, x=1/phi
    mu <- prob %*% t(m)  
    xnew <- x - phifun(x, ctable, mu)/dphifun(x, ctable, mu)
    iternum_phi <- 1
    while (iternum_phi <= iter_max) { 
      x <- xnew
      phi_fun <- phifun(x, ctable, mu)
      if( abs(phi_fun) < 1e-8 ) {break} 
      iternum_phi <- iternum_phi + 1
      xnew <- x - phi_fun / dphifun(x, ctable, mu)
      if ( abs(xnew/x - 1) < 1e-6 ) {break}
      if(xnew > 1e10) { xnew <- Inf; break }
      if(xnew < -1e10) { xnew <- -Inf; break } 
    }
    phinew <- 1/xnew
    if(phinew == 0) { break } 
    if(abs(phinew/phi-1) < 1e-6) { break }
  }
  phi <- phinew
  
  if(phi < 0) {
    phi <- 0
    w <- m/sum(m)
    prob <- rowSums(ctable) / sum(m)
  }
  
  obj <- list(phi, w, as.vector(prob), iternum)
  names(obj) <- c("phi", "weights", "prob", "iternum")
  obj
}

calVar <- function(q, m) {
  p <- q$prob
  phi <- q$phi
  w <- q$weights
  w2 <- w * w
  w2m <- w2 / m
  sum1 <- sum(w2)
  sum2 <- sum(w2m)
  phi * sum1 * p * p + sum2 * p
} 
