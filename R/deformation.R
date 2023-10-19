# .makeE
# .makeXSk
# .makeXS
# .betaxy
# .new_xy
# .reml0
# .d3_deform
# .reml1 analytical
# .reml1 finite-difference
# .d0_deform
# .d12_deform
# .vcov_deformation
# .search_deform
# .makeJ
# .inits
# .predict.deformation
# .plot.deformation
# .beta2delta

.makeE <- function(xy, xy2 = xy) {
  Dsq <- outer(xy2[, 1], xy[, 1], FUN = '-')^2 + outer(xy2[, 2], xy[, 2], FUN = '-')^2
  E <- Dsq * log(sqrt(Dsq))
  E[!is.finite(E)] <- 0
  E
}

.makeXSk <- function(eE, T, k, M) {
  Uk <- eE$vectors[, seq_len(k), drop = FALSE]
  dk <- eE$values[seq_len(k)]
  Zk <- qr.Q(qr(crossprod(Uk, T)), complete = TRUE)
  Zk <- Zk[, ncol(Zk) - (k - M - 1):0, drop = FALSE]
  Xk <- crossprod(t(Uk) * dk, Zk)
  Xk <- cbind(T[, -1], Xk)
  Sk <- diag(numeric(k - 1))
  Sk[2 + 1:(k - M), 2 + 1:(k - M)] <- crossprod(Zk * dk, Zk)
  list(X = Xk, S = Sk, UZ = Uk %*% Zk)
}

.makeXS <- function(eE, T, k, M) {
XS <- lapply(k, .makeXSk, eE = eE, T = T, M = M)
S <- lapply(XS, '[[', 'S')
S <- lapply(S, function(x) x / norm(x, 'F'))
attr(S, 'logdet') <- sapply(S, .gdet)
attr(S, 'rank') <- sapply(S, .rank)
list(X = lapply(XS, '[[', 'X'), S = S, UZ = lapply(XS, '[[', 'UZ'))
}

.betaxy <- function(pars, id) {
  betas <- split(pars, id)
  betas[[1]][1:2] <- exp(betas[[1]][1:2])
  betax <- c(betas[[1]][c(1, 3)], betas[[2]])
  betay <- c(betas[[1]][c(3, 2)], betas[[3]])
  list(betax = betax, betay = betay, beta = betas)
}

.new_xy <- function(pars, id, X, atts = TRUE){
  betaxy <- .betaxy(pars, id)
  out <- cbind(X[[1]] %*% betaxy[[1]], X[[2]] %*% betaxy[[2]])
  if (atts)
    attr(out, 'beta') <- betaxy
  out
}

.reml0 <- function(pars, id, XS, V, m, covfn, lambda = NULL, biject) {
  beta <- attr(pars, 'beta')
  if (is.null(lambda)) {
    lambda <- exp(as.vector(pars))
  } else {
    lambda[attr(lambda, 'fit')] <- exp(as.vector(pars))
  }
  fit <- .newton_step_inner(beta, .d0_deform, .search_deform, id = id, XS = XS, V = V, lambda = lambda, m = m, covfn = covfn, biject = biject, itlim = 1e4)
  out1 <- as.vector(fit$objective)
  out2 <- - .5 * sum(log(lambda) * attr(XS$S, 'rank'))
  out3 <- .5 * .logdetH(fit$hessian[-seq_len(covfn$n0), -seq_len(covfn$n0)])$ldetH
  out <- c(out1, out2, out3)
  out <- sum(out)
  attr(out, 'beta') <- fit$par
  attr(out, 'beta2') <- attr(fit$objective, 'beta')
  attr(out, 'gradient') <- fit$gradient
  attr(out, 'hessian') <- fit$hessian
  out
}

.d3_deform <- function(pars, id, XS, V, lambda, m, covfn, eps = 1e-5) {
H0 <- as.vector(.d12_deform(pars, id, XS, V, lambda, m, covfn)$H)
out <- matrix(0, length(H0), length(pars))
for (i in seq_along(pars)) {
  parsi <- replace(pars, i, pars[i] + eps)
  out[, i] <- .d12_deform(parsi, id, XS, V, lambda, m, covfn)$H
}
(out - H0) / eps
}

.reml1 <- function(pars, id, XS, V, m, covfn, lambda = NULL, biject, eps = 1e-3) {
f0 <- .reml0(pars, id, XS, V, m, covfn, lambda, biject)
b0 <- attr(f0, 'beta')
mult <- c(1, -1)[as.integer(pars > 0) + 1]
oute <- numeric(length(pars))
for (i in seq_along(pars)) {
  pe <- replace(pars, i, pars[i] + mult[i] * eps)
  attr(pe, 'beta') <- b0
  oute[i] <- mult[i] * (.reml0(pe, id, XS, V, m, covfn, lambda, biject) - f0)
}
as.vector(oute) / eps
}

.d0_deform <- function(pars, id, XS, V, lambda, m, covfn, biject = FALSE) {
  pars0 <- pars
  not_beta <- covfn$not_beta
  beta <- .betaxy(pars[-not_beta], id)
  id2 <- c(not_beta, length(not_beta) + rep(seq_len(2), sapply(beta[1:2], length)))
  pars <- c(pars[not_beta], unlist(beta[1:2]))
  C <- covfn$d0(split(pars, id2), XS$X, nrow(V))
  cholC <- try(chol(C), silent = TRUE)
  if (inherits(cholC, 'try-error'))
    return(1e8)
  out1 <- (m - 1) * sum(log(diag(cholC)))
  out2 <- m * sum(diag(backsolve(cholC, backsolve(cholC, V, transpose = TRUE))))
  out <- out1 + .5 * out2
  if (biject) {
    bij_pen <- .d0_biject(unlist(beta[1:2]), id, XS, biject)
    out <- out + bij_pen
  }
  pen <- .5 * sum(lambda * mapply(function(x, y) sum(x * (y %*% x)), beta[1:2], XS$S))
  out <- out + pen
  if (!is.finite(out))
    out <- 1e20
  jeff <- pars[1]
  out <- out + jeff
  attr(out, 'beta') <- beta
  out
}

.d12_deform <- function(pars, id, XS, V, lambda, m, covfn, biject = FALSE) {
  pars0 <- pars
  not_beta <- covfn$not_beta
  beta <- .betaxy(pars[-not_beta], id)
  kk <- sapply(beta[1:2], length)
  J <- .makeJ(length(pars), kk + 1, exp(pars[-not_beta][1:2]), covfn$n0)
  id2 <- c(not_beta, length(not_beta) + rep(seq_len(2), kk))
  pars <- c(pars[not_beta], unlist(beta[1:2]))
  np <- length(id2)
  gH <- .dSigma(pars, XS$X, id, id2, XS$D, covfn, nderiv = 2)
  iAdA <- lapply(gH$d1, function(x) .solve_chol(gH$cholA, x))
  g1 <- (m - 1) * sapply(iAdA, function(x) sum(diag(x)))
  g2 <- m * sapply(gH$d1, function(x) - sum(.solve_chol(gH$cholA, V) * t(.solve_chol(gH$cholA, x))))
  g <- .5 * as.vector(g1 + g2)
  g[-not_beta] <- g[-not_beta] + unlist(Map('*', lambda, Map('%*%', XS$S, beta[1:2])))
  g[1] <- g[1] + 1
  H <- diag(rep(0, np))
  for (i in 1:np) for (j in i:np) {
    temph <- .solve_chol(gH$cholA, gH$d2[[i]][[j]])
    H1 <- (m - 1) * (sum(diag(temph)) - sum(t(iAdA[[i]]) * iAdA[[j]]))
    H2 <- m * sum(diag(.solve_chol(gH$cholA, crossprod(V, 2 * iAdA[[i]] %*% iAdA[[j]] - temph))))
    H[j, i] <- H[i, j] <- .5 * (H1 + H2)
  }
  H <- H + .blockdiag(.clist(list(diag(numeric(length(not_beta))), Map('*', lambda, XS$S))))
  if (biject) {
    gH_biject <- .d12_biject(unlist(beta[1:2]), id, XS, biject)
    g[-not_beta] <- g[-not_beta] + gH_biject$g
    H[-not_beta, -not_beta] <- H[-not_beta, -not_beta] + gH_biject$H
  }
  g <- J %*% g
  H <- tcrossprod(J %*% H, J)
  list(g = g, H = H)
}

.vcov_deformation <- function(pars, id, XS, V, lambda, m, covfn, biject = FALSE) {
  pars0 <- pars
  not_beta <- covfn$not_beta
  beta <- .betaxy(pars[-not_beta], id)
  kk <- sapply(beta[1:2], length)
  # J <- .makeJ(length(pars), kk + 1, exp(pars[5:6]))
  J <- .makeJ(length(pars), kk + 1, exp(pars[-not_beta][1:2]), covfn$n0)
  id2 <- c(not_beta, length(not_beta) + rep(seq_len(2), kk))
  pars <- c(pars[not_beta], unlist(beta[1:2]))
  np <- length(id2)
  gH <- .dSigma(pars, XS$X, id, id2, XS$D, covfn, nderiv = 2)
  iAdA <- lapply(gH$d1, function(x) .solve_chol(gH$cholA, x))
  H <- diag(rep(0, np))
  for (i in 1:np) for (j in i:np) {
    temph <- .solve_chol(gH$cholA, gH$d2[[i]][[j]])
    H1 <- (m - 1) * (sum(diag(temph)) - sum(t(iAdA[[i]]) * iAdA[[j]]))
    H2 <- m * sum(diag(.solve_chol(gH$cholA, crossprod(V, 2 * iAdA[[i]] %*% iAdA[[j]] - temph))))
    H[j, i] <- H[i, j] <- .5 * (H1 + H2)
  }
  H <- H + .blockdiag(.clist(list(diag(numeric(length(not_beta))), Map('*', lambda, XS$S))))
  if (biject) {
    gH_biject <- .d12_biject(unlist(beta[1:2]), id, XS, biject)
  #  g[-not_beta] <- g[-not_beta] + gH_biject$g
    H[-not_beta, -not_beta] <- H[-not_beta, -not_beta] + gH_biject$H
  }
  H <- tcrossprod(J %*% H, J)
  list(H = H, J = J, id = id, Vp = MASS::ginv(H))
}

.search_deform <- function(pars, id, XS, V, lambda, m, covfn, biject) {
gH <- .d12_deform(pars, id, XS, V, lambda, m, covfn, biject)
out <- .solve_chol2(gH$H, gH$g)
attr(out, 'gradient') <- gH[[1]]
attr(out, 'Hessian') <- gH[[2]]
out
}

.makeJ <- function(m, k, eb, n0) {
  ind0 <- c(1:7, 7:(sum(k) + 1))
  ind1 <- c(1:5, k[1] + 5, 6, k[1] + 4, 7:(k[1] + 3), (k[1] + 6):(sum(k) + 2))
  r1 <- c(1:(n0 + 3), n0 + 2 + 1:(sum(k - 3) + 1))
  t2 <- seq_len(k[1] - 3) + n0 + 3
  t3 <- seq_len(k[2] - 3) + k[1] + n0
  r1 <- c(1:n0, n0 + c(1, 3), t2, n0 + c(3, 2), t3)
  c1 <- seq_len(m + 1)
  J <- J2 <- matrix(0, m, m + 1)
  J2[cbind(r1, c1)] <- 1
  J[cbind(5:6, c(5, k[1] + 5))] <- J[cbind(5:6, c(5, k[1] + 5))] * eb
  r2 <- n0 + 1:2
  c2 <- n0 + c(1, k[1] + 1)
  J2[cbind(r2, c2)] <- J2[cbind(r2, c2)] * eb
  J2
}

.inits <- function(x, y, V, m, covfn, hessian = FALSE) {
X <- list(matrix(x, ncol = 1), matrix(y, ncol = 1))
rX <- mean(sapply(X, function(x) diff(range(x))))
inits <- c(.5, log(max(V)), 2, 0, 0)
names(inits) <- c('p0', 'p1', 'p2', 'p3', 'phi1')
inits <- inits[covfn$nms]
inits <- optim(inits, .rss_iso, XS = list(X = X), V = V, m = m, covfn = covfn, control = list(maxit = 1e3))$par
if ('p3' %in% names(inits))
  inits['p3'] <- min(1, inits['p3'])
if ('p2' %in% names(inits))
  inits['p2'] <- min(2, inits['p2'])
inits <- nlm(.d0_iso, inits, XS = list(X = X), V = V, m = m, covfn = covfn)$estimate
inits <- c(inits, inits[length(inits)])
opt <- optim(inits, .d0_aniso, XS = list(X = X), V = V, m = m, covfn = covfn, control = list(maxit = 1e3), hessian = hessian)
out <- opt$par
if ('p3' %in% names(out))
  out['p3'] <- min(1, out['p3'])
if ('p2' %in% names(inits))
  out['p2'] <- min(2, out['p2'])
if ('p0' %in% names(inits))
  out['p0'] <- min(3, out['p0'])
if (hessian) {
  H <- opt$hessian
  attr(out, 'hessian') <- H
}
out
}

.std_fn <- function(x, scl = NULL, common = TRUE) {
  if (is.null(scl)) {
    scl <- colMeans(x)
    if (common) {
      scl <- rbind(scl, mean(apply(x, 2, sd)))
    } else {
      scl <- rbind(scl, apply(x, 2, sd))
    }  
  }
  out <- t((t(x) - scl[1, ]) / scl[2, ])
  attr(out, 'scaling') <- scl
  out
}

.predict.deformation <- function(object, newx, se.fit = FALSE) {
if (is.null(newx)) {
  out <- object$fitted
  E <- object$E
  newx <- object$xy
} else {
  if (!is.null(object$scaling)) {
    newx <- .std_fn(newx, object$scaling)
  }
  E <- .makeE(object$xy, newx)
  T <- cbind(1, newx)
  out <- newx %*% cbind(object$beta[[1]][1:2], object$beta[[2]][1:2])
  out <- out + E %*% object$delta
}
if (se.fit) {
  out <- list(fitted = out)
  Vp <- crossprod(object$V$J, object$V$Vp %*% object$V$J)[-object$not_beta, -object$not_beta]
  X <- lapply(object$XS$UZ, function(x) cbind(newx, E %*% x))
  id <- rep(seq_along(X), sapply(X, ncol))
  Vp <- lapply(seq_along(X), function(i) Vp[id == i, id == i])
  ese <- mapply(function(x, y) rowSums(x * t(tcrossprod(y, x))), X, Vp)
  out$se.fit <- sqrt(ese)
}
out
}

.plot.deformation <- function(x, nx, ny, xp, yp, xlab, ylab, ...) {
if (is.null(xp))
  xp <- pretty(x$xy0[, 1], nx)
if (is.null(yp))
  yp <- pretty(x$xy0[, 2], ny)
if (is.null(xlab))
  xlab <- 'x*'
if (is.null(ylab))
  ylab <- 'y*'
xyp <- as.matrix(expand.grid(xp, yp))
xyp2 <- .predict.deformation(x, xyp)
xm <- matrix(xyp2[, 1], length(xp))
ym <- matrix(xyp2[, 2], length(xp))
matplot(xm, ym, type = 'l', col = 1, lty = 1, xlab = xlab, ylab = ylab, ...)
matlines(t(xm), t(ym), col = 1, lty = 1)
}

.beta2delta <- function(beta, UZl) {
sapply(1:2, function(i) as.vector(UZl[[i]] %*% beta[[3]][[i + 1]]))
}
