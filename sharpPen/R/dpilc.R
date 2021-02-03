dpilc <- function(xx, yy, blockmax = 5, divisor = 20, trim = 0.01,
                  proptrun = 0.05, gridsize = 401L,
                  range.x = range(x), truncate = TRUE)
{
  ## Trim the 100(trim)% of the data from each end (in the x-direction).
  blkest46 <- function(x, y, Nval, q)
  {
    n <- length(x)
    
    ## Sort the (x, y) data with respect to
    ## the x's.
    
    datmat <- cbind(x, y)
    datmat <- datmat[sort.list(datmat[, 1L]), ]
    x <- datmat[, 1L]
    y <- datmat[, 2L]
    
    ## Set up arrays for FORTRAN programme "blkest"
    
    qq <- q + 1L
    xj <- rep(0, n)
    yj <- rep(0, n)
    coef <- rep(0, qq)
    Xmat <- matrix(0, n, qq)
    wk <- rep(0, n)
    qraux <- rep(0, qq)
    sigsqe <- 0
    th44e <- 0
    th46e <- 0
    
    out <- .Fortran("blkest46", as.double(x), as.double(y), as.integer(n),
                    as.integer(q), as.integer(qq), as.integer(Nval), as.double(xj),
                    as.double(yj), as.double(coef), as.double(Xmat), as.double(wk),
                    as.double(qraux), as.double(sigsqe), as.double(th44e),
                    as.double(th46e), PACKAGE = "sharpPen")
    
    list(sigsqe = out[[13]], th44e = out[[14]], th46e = out[[15]])
  }
  
  
  cpblock <- function(X, Y, Nmax, q)
  {
    n <- length(X)
    
    ## Sort the (X, Y) data with respect tothe X's.
    
    datmat <- cbind(X, Y)
    datmat <- datmat[sort.list(datmat[, 1L]), ]
    X <- datmat[, 1L]
    Y <- datmat[, 2L]
    
    ## Set up arrays for FORTRAN subroutine "cp"
    
    qq <- q + 1L
    RSS <- rep(0, Nmax)
    Xj <- rep(0, n)
    Yj <- rep(0, n)
    coef <- rep(0, qq)
    Xmat <- matrix(0, n, qq)
    Cpvals <- rep(0, Nmax)
    wk <- rep(0, n)
    qraux <- rep(0, qq)
    
    ## remove unused 'q' 2007-07-10
    out <- .Fortran("cpsd", as.double(X), as.double(Y), as.integer(n),
                    as.integer(qq), as.integer(Nmax), as.double(RSS), as.double(Xj),
                    as.double(Yj), as.double(coef), as.double(Xmat), as.double(wk),
                    as.double(qraux), Cpvals = as.double(Cpvals), PACKAGE = "sharpPen")
    
    Cpvec <- out$Cpvals
    
    order(Cpvec)[1L]
  }
  
  
  
  rlbin <- function(X, Y, gpoints, truncate = TRUE)
  {
    n <- length(X)
    M <- length(gpoints)
    trun <- if (truncate) 1L else 0L
    a <- gpoints[1L]
    b <- gpoints[M]
    out <- .Fortran("rlbinsd", as.double(X), as.double(Y), as.integer(n),
                    as.double(a), as.double(b), as.integer(M), as.integer(trun),
                    double(M), double(M), PACKAGE = "sharpPen")
    list(xcounts = out[[8L]],  ycounts = out[[9L]])
  }
  
  
  linbin <- function(X, gpoints, truncate = TRUE)
  {
    n <- length(X)
    M <- length(gpoints)
    trun <- if (truncate) 1L else 0L
    a <- gpoints[1L]
    b <- gpoints[M]
    .Fortran("linbinsd", as.double(X), as.integer(n),
             as.double(a), as.double(b), as.integer(M),
             as.integer(trun), double(M), PACKAGE = "sharpPen")[[7]]
  }
  
  
  sdiag <- function(x, drv = 0L, degree = 1L, kernel = "normal",
                    bandwidth, gridsize = 401L, bwdisc = 25, range.x,
                    binned = FALSE, truncate = TRUE)
  {
    if (missing(range.x) && !binned) range.x <- c(min(x), max(x))
    
    ## Rename common variables
    
    M <- gridsize
    Q <- as.integer(bwdisc)
    a <- range.x[1L]
    b <- range.x[2L]
    pp <- degree + 1L
    ppp <- 2L*degree + 1L
    tau <- 4
    
    ## Bin the data if not already binned
    
    if (!binned) {
      gpoints <- seq(a, b, length = M)
      xcounts <- linbin(x, gpoints, truncate)
    } else {
      xcounts <- x
      M <- length(xcounts)
      gpoints <- seq(a, b, length = M)
    }
    
    ## Set the bin width
    
    delta <- (b-a)/(M-1L)
    
    ## Discretise the bandwidths
    
    if (length(bandwidth) == M) {
      hlow <- sort(bandwidth)[1L]
      hupp <- sort(bandwidth)[M]
      hdisc <- exp(seq(log(hlow), log(hupp), length = Q))
      
      ## Determine value of L for each member of "hdisc"
      Lvec <- floor(tau*hdisc/delta)
      
      ## Determine index of closest entry of "hdisc"
      ## to each member of "bandwidth"
      indic <- if (Q > 1L) {
        lhdisc <- log(hdisc)
        gap <- (lhdisc[Q]-lhdisc[1L])/(Q-1)
        if (gap == 0) rep(1, M)
        else round(((log(bandwidth) - log(sort(bandwidth)[1L]))/gap) + 1)
      } else rep(1, M)
    } else if (length(bandwidth) == 1L) {
      indic <- rep(1, M)
      Q <- 1L
      Lvec <- rep(floor(tau*bandwidth/delta), Q)
      hdisc <- rep(bandwidth, Q)
    } else
      stop("'bandwidth' must be a scalar or an array of length 'gridsize'")
    
    dimfkap <- 2L * sum(Lvec) + Q
    fkap <- rep(0, dimfkap)
    midpts <- rep(0, Q)
    ss <- matrix(0, M, ppp)
    Smat <- matrix(0, pp, pp)
    work <- rep(0, pp)
    det <- rep(0, 2L)
    ipvt <- rep(0, pp)
    Sdg <- rep(0, M)
    
    out <- .Fortran("sdiagsd", as.double(xcounts), as.double(delta),
                    as.double(hdisc), as.integer(Lvec), as.integer(indic),
                    as.integer(midpts), as.integer(M), as.integer(Q),
                    as.double(fkap), as.integer(pp), as.integer(ppp),
                    as.double(ss), as.double(Smat), as.double(work),
                    as.double(det), as.integer(ipvt), as.double(Sdg), PACKAGE = "sharpPen")
    
    list(x = gpoints,  y = out[[17L]])
  }
  
  
  
  sstdiag <- function(x, drv = 0L, degree = 1L, kernel = "normal",
                      bandwidth, gridsize = 401L, bwdisc = 25, range.x,
                      binned = FALSE, truncate = TRUE)
  {
    if (missing(range.x) && !binned) range.x <- c(min(x), max(x))
    
    ## Rename common variables
    M <- gridsize
    Q <- as.integer(bwdisc)
    a <- range.x[1L]
    b <- range.x[2L]
    pp <- degree + 1L
    ppp <- 2L*degree + 1L
    tau <- 4L
    
    ## Bin the data if not already binned
    if (!binned) {
      gpoints <- seq(a, b, length = M)
      xcounts <- linbin(x, gpoints, truncate)
    } else {
      xcounts <- x
      M <- length(xcounts)
      gpoints <- seq(a, b, length = M)
    }
    
    ## Set the bin width
    
    delta <- (b-a)/(M-1L)
    
    ## Discretise the bandwidths
    if (length(bandwidth) == M) {
      hlow <- sort(bandwidth)[1L]
      hupp <- sort(bandwidth)[M]
      hdisc <- exp(seq(log(hlow), log(hupp), length = Q))
      
      ## Determine value of L for each member of "hdisc"
      Lvec <- floor(tau*hdisc/delta)
      
      ## Determine index of closest entry of "hdisc"
      ## to each member of "bandwidth"
      indic <- if (Q > 1L) {
        lhdisc <- log(hdisc)
        gap <- (lhdisc[Q]-lhdisc[1L])/(Q-1)
        if (gap == 0) rep(1, M)
        else round(((log(bandwidth) - log(sort(bandwidth)[1L]))/gap) + 1)
      } else rep(1, M)
    } else if (length(bandwidth) == 1L) {
      indic <- rep(1, M)
      Q <- 1L
      Lvec <- rep(floor(tau*bandwidth/delta), Q)
      hdisc <- rep(bandwidth, Q)
    } else
      stop("'bandwidth' must be a scalar or an array of length 'gridsize'")
    
    dimfkap <- 2L * sum(Lvec) + Q
    fkap <- rep(0, dimfkap)
    midpts <- rep(0, Q)
    ss <- matrix(0, M, ppp)
    uu <- matrix(0, M, ppp)
    Smat <- matrix(0, pp, pp)
    Umat <- matrix(0, pp, pp)
    work <- rep(0, pp)
    det <- rep(0, 2L)
    ipvt <- rep(0, pp)
    SSTd <- rep(0, M)
    
    SSTd <- .Fortran("sstdgsd", as.double(xcounts), as.double(delta),
                     as.double(hdisc), as.integer(Lvec), as.integer(indic),
                     as.integer(midpts), as.integer(M), as.integer(Q),
                     as.double(fkap), as.integer(pp), as.integer(ppp),
                     as.double(ss), as.double(uu), as.double(Smat),
                     as.double(Umat), as.double(work), as.double(det),
                     as.integer(ipvt), as.double(SSTd), PACKAGE = "sharpPen")[[19L]]
    
    list(x = gpoints, y = SSTd)
  }
  
  
  xy <- cbind(xx, yy)
  xy <- xy[sort.list(xy[, 1L]), ]
  x <- xy[, 1L]
  y <- xy[, 2L]
  indlow <- floor(trim*length(x)) + 1
  indupp <- length(x) - floor(trim*length(x))
  
  x <- x[indlow:indupp]
  y <- y[indlow:indupp]
  
  ## Rename common parameters
  n <- length(x)
  M <- gridsize
  a <- range.x[1L]
  b <- range.x[2L]
  
  ## Bin the data
  
  gpoints <- seq(a, b, length = M)
  out <- rlbin(x, y, gpoints, truncate)
  xcounts <- out$xcounts
  ycounts <- out$ycounts
  
  ## Choose the value of N using Mallow's C_p
  Nmax <- max(min(floor(n/divisor), blockmax), 1)
  Nval <- cpblock(x, y, Nmax, 6)
  
  ## Estimate sig^2, theta_22 and theta_24 using quartic fits
  ## on "Nval" blocks.
  
  out <- blkest46(x, y, Nval, 6)
  sigsqQ <- out$sigsqe
  th46Q <- out$th46e

  ## Estimate theta_22 using a local cubic fit
  ## with a "rule-of-thumb" bandwidth: "gamseh"
  
  gamseh <- (sigsqQ*(b-a)/(abs(th46Q)*n))
  Rk45<-105/(32*sqrt(pi))
  uk45<-360
  if (th46Q < 0) gamseh <- (((44-28)*22.5)*(Rk45/uk45)*gamseh)^(1/11)
  if (th46Q > 0) gamseh <- (((44+28)*22.5)*(Rk45/uk45)*gamseh)^(1/11)
  
  mddest <- locpoly(xcounts, ycounts, drv=4L, bandwidth=gamseh,
                    range.x=range.x, binned=TRUE)$y
  
  llow <- floor(proptrun*M) + 1
  lupp <- M - floor(proptrun*M)
  th44kn <- sum((mddest[llow:lupp]^2)*xcounts[llow:lupp])/n
  ## Estimate sigma^2 using a local linear fit
  ## with a "direct plug-in" bandwidth: "lamseh"
  C3K0 <- 7881/(32*16^2*sqrt(2*pi))+27/(8*sqrt(pi))-4192/(9*64*sqrt(6*pi))
  C3K <- (8^3*C3K0)^(1/17)
  lamseh <- C3K*(((sigsqQ^2)*(b-a)/((th44kn*n)^2))^(1/17))
  
  ## Now compute a local linear kernel estimate of
  ## the variance.
  mest <- locpoly(xcounts, ycounts, bandwidth=lamseh,
                  range.x=range.x, binned=TRUE)$y
  Sdg <- sdiag(xcounts, bandwidth=lamseh,
               range.x=range.x, binned=TRUE)$y
  SSTdg <- sstdiag(xcounts, bandwidth=lamseh,
                   range.x=range.x, binned=TRUE)$y
  sigsqn <- sum(y^2) - 2*sum(mest*ycounts) + sum((mest^2)*xcounts)
  sigsqd <- n - 2*sum(Sdg*xcounts) + sum(SSTdg*xcounts)
  sigsqkn <- sigsqn/sigsqd
  
  ## Combine to obtain final answer.
  result<-(27*sigsqkn*(b-a)/(4*sqrt(pi)*th44kn*n))^(1/9)
  result
  
}
