# .d0_iso
# .d0_aniso
# .predict.aniso
# .plot.aniso

.d0_iso <- function(pars, id, XS, V, m, covfn) {
  if (pars[3] > 5)
    return(1e8)
  not_beta <- covfn$not_beta
  pars <- c(pars[not_beta], rep(exp(pars[-not_beta]), 2))
  id2 <- seq_along(pars)
  C <- covfn$d0(split(pars, id2), XS$X, nrow(V))
  cholC2 <- try(chol(C, pivot = TRUE), silent = TRUE)
  if (attr(cholC2, 'rank') < ncol(V))
    return(1e8)
  piv <- attr(cholC2, 'pivot')
  ipiv <- order(piv)
  out1 <- (m - 1) * sum(log(diag(cholC2)))
  out2 <- m * sum(diag(backsolve(cholC2, forwardsolve(cholC2, V[piv, piv], transpose = TRUE, upper.tri = TRUE))[ipiv, ipiv]))
  out <- out1 + .5 * out2
  if (!is.finite(out))
    out <- 1e20
  out
}

.d0_aniso <- function(pars, id, XS, V, m, covfn) {
  if (pars[3] > 5)
    return(1e8)
  not_beta <- covfn$not_beta
  pars <- c(pars[not_beta], exp(pars[-not_beta]))
  id2 <- seq_along(pars)
  C <- covfn$d0(split(pars, id2), XS$X, nrow(V))
  cholC2 <- try(chol(C, pivot = TRUE), silent = TRUE)
  if (attr(cholC2, 'rank') < ncol(V))
    return(1e8)
  piv <- attr(cholC2, 'pivot')
  ipiv <- order(piv)
  out1 <- (m - 1) * sum(log(diag(cholC2)))
  out2 <- m * sum(diag(backsolve(cholC2, forwardsolve(cholC2, V[piv, piv], transpose = TRUE, upper.tri = TRUE))[ipiv, ipiv]))
  out <- out1 + .5 * out2
  if (!is.finite(out))
    out <- 1e20
  out
}

.predict.aniso <- function(object, newx, se.fit = FALSE) {
if (is.null(newx)) {
  out <- object$fitted
} else {
  if (!is.null(object$scaling)) {
    newx <- .std_fn(newx, object$scaling)
  }
  out <- t(object$b * t(newx))
}
if (se.fit) {
  out <- list(fitted = out)
  ese <- sqrt(object$b * tail(diag(object$V$Vp), 2))
  out$se.fit <- matrix(ese, nrow(out$fitted), 2, byrow = TRUE)
}
out
}

.plot.aniso <- function(x, nx, ny, xp, yp, xlab, ylab, ...) {
if (is.null(xp))
  xp <- pretty(x$xy[, 1], nx)
if (is.null(yp))
  yp <- pretty(x$xy[, 2], ny)
if (is.null(xlab))
  xlab <- 'x*'
if (is.null(ylab))
  ylab <- 'y*'
xyp <- as.matrix(expand.grid(xp, yp))
xyp2 <- .predict.aniso(x, xyp)
xm <- matrix(xyp2[, 1], length(xp))
ym <- matrix(xyp2[, 2], length(xp))
matplot(xm, ym, type = 'l', col = 1, lty = 1, xlab = xlab, ylab = ylab, ...)
matlines(t(xm), t(ym), col = 1, lty = 1)
}

.rss_iso <- function(pars, id, XS, V, m, covfn) {
  if (pars[3] > 5)
    return(1e8)
  not_beta <- covfn$not_beta
  pars <- c(pars[not_beta], rep(exp(pars[-not_beta]), 2))
  id2 <- seq_along(pars)
  C <- covfn$d0(split(pars, id2), XS$X, nrow(V))
  rss <- C - V
  sum(rss[!lower.tri(rss)]^2)
}

.rss_aniso <- function(pars, id, XS, V, m, covfn) {
  if (pars[3] > 5)
    return(1e8)
  not_beta <- covfn$not_beta
  pars <- c(pars[not_beta], exp(pars[-not_beta]))
  id2 <- seq_along(pars)
  C <- covfn$d0(split(pars, id2), XS$X, nrow(V))
  rss <- C - V
  sum(rss[!lower.tri(rss)]^2)
}
