# .gdet
# .rank
# .reml0
# .longlist
# .clist
# .vecdiff
# .dSigma
# .solve_chol2
# .doPD
# .semivariog
# .inits
# .blockdiag
# .make_cov_fns
# .powexp_fns
# .unitpowexp_fns
# .cospowexp_fns
# .unitcospowexp_fns

.gdet <- function(x, tol = .Machine$double.eps^.7, log = TRUE) {
cd <- suppressWarnings(chol(x, pivot = TRUE))
dx <- diag(cd)
out <- 2 * sum(log(dx[dx > tol]))
if (!log)
  out <- exp(out)
out
}

.rank <- function(x, tol = .Machine$double.eps^.7, log = TRUE) {
cd <- suppressWarnings(chol(x, pivot = TRUE))
attr(cd, 'rank')
}

.longlist <- function(x) {
out <- list()
it <- 1
for (i in seq_along(x)) {
if (inherits(x[[i]], "list")) {
for (j in seq_along(x[[i]])) {
out[[it]] <- x[[i]][[j]]
it <- it + 1
}
} else {
out[[it]] <- x[[i]]
it <- it + 1
}
}
out
}

.clist <- function(x) {
  x <- .longlist(x)
  while(any(sapply(x, is.list))) {
    x <- .longlist(x)
  }
  x
}

.vecdiff <- function(x) t(outer(x, x, FUN = '-'))

.dSigma <- function(pars, X, id, id2, DXxy, covfn, nderiv = 0) {

ndeform <- length(X)
ncov <- max(id2) - ndeform
np <- length(id2)

A0 <- covfn$d0(split(pars, id2), X, nrow(X[[1]]))

if (any(!is.finite(A0)))
  stop("Some non-finite covariance matrix entries.")

out <- list(d0 = A0)

if (nderiv >= 1) {

idX1 <- rep(seq_along(DXxy), sapply(DXxy, length))
idX2 <- unlist(lapply(sapply(DXxy, length), seq_len))

out$cholA <- suppressWarnings(chol(A0, pivot = TRUE))
out$i0 <- .solve_chol(out$cholA, diag(nrow(A0)))

dA <- covfn$d1(split(pars, id2), X, nrow(X[[1]]))

for (i in ncov + seq_len(ndeform))
  dA[[i]] <- lapply(DXxy[[i - ncov]], function(x) x * dA[[i]])

dA <- .clist(dA)

out$d1 <- dA

if (nderiv >= 2) {

d2A <- covfn$d2(split(pars, id2), X, nrow(X[[1]]))
d2A2 <- lapply(seq_along(pars), function(i) lapply(seq_along(pars), function(x) NULL))

for (i in 1:np) for (j in i:np) {
  temph <- d2A[[id2[i]]][[id2[j]]]
  if (id2[j] > ncov) {
    idXj <- c(idX1[[j - ncov]], idX2[[j - ncov]]) # position in DXxy
    if (id2[i] <= ncov) {
      temph <- temph * DXxy[[idXj[1]]][[idXj[2]]]
    } else {
      idXi <- c(idX1[[i - ncov]], idX2[[i - ncov]]) # position in DXxy
      temph <- temph * DXxy[[idXi[1]]][[idXi[2]]] * DXxy[[idXj[1]]][[idXj[2]]]
    }
  }
  temph[!is.finite(temph)] <- 0
  d2A2[[i]][[j]] <- temph
}

out$d2 <- d2A2

}

}

out

}

.solve_chol2 <- function(H, x) {
d <- 1 / sqrt(diag(abs(H)))
D <- diag(d)
H2 <- D %*% H %*% D
A <- A0 <- H2
eps <- 1e-12
cA <- suppressWarnings(chol(A, pivot = TRUE))
while (attr(cA, 'rank') < nrow(A)) {
  A <- A0
  diag(A) <- diag(A) + eps
  cA <- suppressWarnings(chol(A, pivot = TRUE))
  eps <- 10 * eps
}
L <- cA
x <- d * x
piv <- ipiv <- attr(L, "pivot")
if (is.null(piv)) {
  out <- backsolve(L, backsolve(L, x, transpose = TRUE))
} else {
  ipiv[piv] <- seq_along(piv)
  out <- d * (backsolve(L, backsolve(L, x[piv, , drop=FALSE], upper.tri=TRUE, transpose=TRUE))[ipiv, , drop=FALSE])
}
attr(out, 'cholH') <- L
attr(out, 'rank') <- attr(L, 'rank')
out
}

.doPD <- function(A) {
A0 <- A
eps <- 1e-12
cA <- suppressWarnings(chol(A, pivot = TRUE))
while (attr(cA, 'rank') < nrow(A)) {
  A <- A0
  diag(A) <- diag(A) + eps
  cA <- suppressWarnings(chol(A, pivot = TRUE))
  eps <- 10 * eps
}
return(cA)
}

.sqdist <- function(x) {
  outer(x, x, FUN = '-')^2
}

.dist2 <- function(x) {
matrix(sqrt(rowSums(apply(x, 2, .sqdist))), nrow(x))
}

.semivariog <- function(x, V, bins = 0, bin.function, trim) {
use <- !upper.tri(V)
d <- .dist2(x)[use]
gamma <- outer(diag(V), diag(V), FUN = '+') - 2 * V
gamma <- gamma[use]
wts <- rep(1, sum(use))
if (bins > 0) {
  if (bin.function == 'quantile') {
    brks <- quantile(d, ppoints(bins))
    brks <- c(-1e-6, brks)
  } else {
    brks <- seq(-1e-6, max(d), l = bins + 1)
  }
  dc <- cut(d, brks)
  spl <- split(gamma, dc)
  gamma <- sapply(spl, mean, trim = trim)
  d <- brks[-1] - .5 * diff(brks)
  wts <- sapply(spl, length)
}
cbind(d, .5 * gamma, wts)
}

.blockdiag <- function(lst) {
ends <- sapply(lst, dim)
starts <- 1 + apply(cbind(0, ends[, seq_len(ncol(ends) - 1)]), 1, cumsum)
ends <- apply(ends, 1, cumsum)
nrc <- ends[length(lst), ]
out <- matrix(0, nrc[1], nrc[2])
for (i in seq_along(lst)) {
  out[starts[i, 1]:ends[i, 1], starts[i, 2]:ends[i, 2]] <- lst[[i]]
}
out
}

.make_cov_fns <- function(correlation, cosine) {
  if (correlation) {
    if (cosine) {
      out <- .unitcospowexp_fns
    } else {
      out <- .unitpowexp_fns
    }
  } else {
    if (cosine) {
      out <- .cospowexp_fns
    } else {
      out <- .powexp_fns
    }
  }
  out
}
  
.powexp_fns <- list(d0 = .covfn, d1 = .d1covfn, d2 = .d2covfn, pd0 = .pcovfn, 
                   nms = c('p1', 'p2', 'p3', 'phi1'), 
                   not_beta = 1:3, n0 = 3)
.unitpowexp_fns <- list(d0 = .unitcovfn, d1 = .d1unitcovfn, d2 = .d2unitcovfn, pd0 = .punitcovfn, 
                       nms = c('p2', 'p3', 'phi1'),
                       not_beta = 1:2, n0 = 2)
.cospowexp_fns <- list(d0 = .coscovfn, d1 = .d1coscovfn, d2 = .d2coscovfn, pd0 = .pcoscovfn, 
                      nms = c('p0', 'p1', 'p2', 'p3', 'phi1'),
                      not_beta = 1:4, n0 = 4)
.unitcospowexp_fns <- list(d0 = .unitcoscovfn, d1 = .d1unitcoscovfn, d2 = .d2unitcoscovfn, pd0 = .punitcoscovfn, 
                          nms = c('p0', 'p2', 'p3', 'phi1'),
                          not_beta = 1:3, n0 = 3)
.dampedcos_fns <- list(d0 = .dampedcoscovfn, d1 = .d1dampedcoscovfn, d2 = .d2dampedcoscovfn, pd0 = .pdampedcoscovfn, 
                       nms = c('p0', 'p1', 'p2', 'phi1'),
                       not_beta = 1:3, n0 = 3)

