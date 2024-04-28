##################################################################
##################################################################
calculateXI <- function (xvec, yvec, simple = TRUE, seed = 12133331) 
{
  #set.seed(seed)  # don't set seed for simulations
  n <- length(xvec)
  PI <- rank(xvec, ties.method = "random")
  fr <- rank(yvec, ties.method = "max")/n
  gr <- rank((-yvec), ties.method = "max")/n
  ord <- order(PI)
  fr <- fr[ord]
  A1 <- sum(abs(fr[1:(n - 1)] - fr[2:n]))/(2 * n)
  CU <- mean(gr * (1 - gr))
  xi <- 1 - A1/CU
  if (simple == TRUE) 
    return(xi)
  else return(list(xi = xi, fr = fr, CU = CU))
}

xicor <- function (x, y = NULL, pvalue = FALSE, ties = TRUE, method = "asymptotic", 
          nperm = 1000, factor = FALSE) 
{
  if (factor == TRUE) {
    if (!is.numeric(x)) 
      x <- as.numeric(factor(x))
    if (!is.numeric(y)) 
      y <- as.numeric(factor(y))
  }
  if (is.data.frame(y)) 
    y <- as.matrix(y)
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.matrix(x) && is.null(y)) 
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y))) 
      stop("'y' must be numeric")
    stopifnot(is.atomic(y))
  }
  if (is.null(y)) {
    ncy <- ncx <- ncol(x)
    if (ncx == 0) 
      stop("'x' is empty")
    if (pvalue == TRUE) 
      stop("testing is not available for matrices")
    r <- matrix(0, nrow = ncx, ncol = ncy)
    for (i in seq_len(ncx)) {
      for (j in seq_len(i)) {
        x2 <- x[, i]
        y2 <- x[, j]
        ok <- complete.cases(x2, y2)
        x2 <- x2[ok]
        y2 <- y2[ok]
        if (any(ok)) {
          r[i, j] <- calculateXI(x2, y2, simple = TRUE)
          r[j, i] <- calculateXI(y2, x2, simple = TRUE)
        }
        else NA
      }
    }
    rownames(r) <- colnames(x)
    colnames(r) <- colnames(x)
    return(r)
  }
  else if (ncol(as.matrix(x)) == 1 & ncol(as.matrix(y)) == 
           1) {
    ok <- complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]
    res <- calculateXI(x, y, simple = FALSE)
    xi <- res$xi
    CU <- res$CU
    n <- length(x)
  }
  else if (ncol(as.matrix(x)) > 1 & ncol(as.matrix(y)) == 1) {
    ok <- complete.cases(cbind(x, y))
    x <- x[ok]
    y <- y[ok]
    res <- calculateXI(x, y, simple = FALSE)
    xi <- res$xi
    CU <- res$CU
    n <- length(x)
  }
  if (pvalue) {
    if (ties == FALSE) 
      return(list(xi = xi, sd = sqrt(2/(5 * n)), pval = 1 - 
                    pnorm(sqrt(n) * xi/sqrt(2/5))))
    if (!(method %in% c("asymptotic", "permutation"))) 
      stop("method for test can only be asymptotic or permutation")
    if (method == "asymptotic") {
      fr <- res$fr
      qfr <- sort(fr)
      ind <- c(1:n)
      ind2 <- 2 * n - 2 * ind + 1
      ai <- mean(ind2 * qfr * qfr)/n
      ci <- mean(ind2 * qfr)/n
      cq <- cumsum(qfr)
      m <- (cq + (n - ind) * qfr)/n
      b <- mean(m^2)
      v <- (ai - 2 * b + ci^2)/(CU^2)
      return(list(xi = xi, sd = sqrt(v/n), pval = 1 - pnorm(sqrt(n) * 
                                                              xi/sqrt(v))))
    }
    if (method == "permutation") {
      rp <- rep(0, nperm)
      for (i in 1:nperm) {
        x1 <- runif(n, 0, 1)
        rp[i] <- calculateXI(x1, y)
      }
      return(list(xi = xi, sd = sqrt(var(rp)), pval = mean(rp > 
                                                             xi)))
    }
  }
  else return(xi)
}
##################################################################
##################################################################
