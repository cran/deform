.makeXSk_expansion <- function(eE, k, xy) {
out <- list()
  if (k == 0) {
    out$X <- xy
  } else {
    Uk <- eE$vectors[, seq_len(k), drop = FALSE]
    out$Uk <- Uk
    dk <- eE$values[seq_len(k)]
    out$X <- t(t(Uk) * dk)
    out$S <- diag(dk)
  }
out
}

.makeXS_expansion <- function(eE, k, xy) {
XS <- lapply(c(0, 0, k), .makeXSk_expansion, eE = eE, xy = xy)
S <- lapply(XS[-c(1:2)], '[[', 'S')
S <- lapply(S, function(x) x / norm(x, 'F'))
Uk <- lapply(XS[-c(1:2)], '[[', 'Uk')
attr(S, 'logdet') <- sapply(S, .gdet)
attr(S, 'rank') <- sapply(S, .rank)
list(X = lapply(XS, '[[', 'X'), S = S, Uk = Uk)
}

.betaz <- function(pars, id) {
  pars[1:2] <- exp(pars[1:2])
  pars <- c(pars[c(1, 3, 3, 2)], pars[-c(1:3)])
  split(pars, id)
}

.reml0_expansion <- function(pars, id, XS, V, m, covfn, lambda = NULL, biject) {
  beta <- attr(pars, 'beta')
  if (is.null(lambda)) {
    lambda <- exp(as.vector(pars))
  } else {
    lambda[attr(lambda, 'fit')] <- exp(as.vector(pars))
  }
  fit <- .newton_step_inner(beta, .d0_expansion, .search_expansion, id = id, XS = XS, V = V, lambda = lambda, m = m, covfn = covfn, itlim = 1e2)
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

.reml1_expansion <- function(pars, id, XS, V, m, covfn = covfn, lambda = NULL, eps = 1e-3) {
f0 <- .reml0_expansion(pars, id, XS, V, m, covfn, lambda)
b0 <- attr(f0, 'beta')
mult <- c(1, -1)[as.integer(pars > 0) + 1]
oute <- numeric(length(pars))
for (i in seq_along(pars)) {
  pe <- replace(pars, i, pars[i] + mult[i] * eps)
  attr(pe, 'beta') <- b0
  oute[i] <- mult[i] * (.reml0_expansion(pe, id, XS, V, m, covfn, lambda) - f0)
}
as.vector(oute) / eps
}

.d0_expansion <- function(pars, id, XS, V, lambda, m, covfn) {
  not_beta <- covfn$not_beta
  beta <- .betaz(pars[-not_beta], id[-not_beta])
  pars <- c(pars[not_beta], unlist(beta))
  C <- covfn$d0(split(pars, id), XS$X, nrow(V))
  cholC <- try(chol(C), silent = TRUE)
  if (inherits(cholC, 'try-error'))
    return(1e8)
  out1 <- (m - 1) * sum(log(diag(cholC)))
  out2 <- m * sum(diag(backsolve(cholC, backsolve(cholC, V, transpose = TRUE))))
  out <- out1 + .5 * out2
  pen <- .5 * sum(lambda * mapply(function(x, y) sum(x * (y %*% x)), beta[-c(1:2)], XS$S))
  out <- out + pen
  if (!is.finite(out))
    out <- 1e20
  jeff <- pars[1]
  out <- out + jeff
  attr(out, 'beta') <- beta
  out
}

.d12_expansion <- function(pars, id, XS, V, lambda, m, covfn) {
  not_beta <- covfn$not_beta
  not_z <- 1:(4 + covfn$n0)
  beta <- .betaz(pars[-not_beta], id[-not_beta])
  J <- rbind(diag(exp(pars[-not_beta][1:2])), 0)
  J <- cbind(J[, 1], c(0, 0, 1), c(0, 0, 1), J[, 2])
  J <- list(diag(length(not_beta)), J, lapply(sapply(beta[-c(1, 2)], length), diag))
  J <- .blockdiag(.clist(J))
  pars <- c(pars[not_beta], unlist(beta))
  np <- length(pars)
  gH <- .dSigma(pars, XS$X, 0, id, XS$D, covfn, nderiv = 2)
  iAdA <- lapply(gH$d1, function(x) .solve_chol(gH$cholA, x))
  g1 <- (m - 1) * sapply(iAdA, function(x) sum(diag(x)))
  g2 <- m * sapply(gH$d1, function(x) - sum(.solve_chol(gH$cholA, V) * t(.solve_chol(gH$cholA, x))))
  g <- .5 * as.vector(g1 + g2)
  g[-not_z] <- g[-not_z] + unlist(Map('*', lambda, Map('%*%', XS$S, beta[-c(1:2)])))
  g[1] <- g[1] + 1
  H <- diag(rep(0, np))
  for (i in 1:np) for (j in i:np) {
    temph <- .solve_chol(gH$cholA, gH$d2[[i]][[j]])
    H1 <- (m - 1) * (sum(diag(temph)) - sum(t(iAdA[[i]]) * iAdA[[j]]))
    H2 <- m * sum(diag(.solve_chol(gH$cholA, crossprod(V, 2 * iAdA[[i]] %*% iAdA[[j]] - temph))))
    H[j, i] <- H[i, j] <- .5 * (H1 + H2)
  }
  H <- H + .blockdiag(.clist(list(diag(numeric(4 + covfn$n0)), Map('*', lambda, XS$S))))
  g <- J %*% g
  H <- tcrossprod(J %*% H, J)
  list(g = g, H = H)
}

.vcov_expansion <- function(pars, id, XS, V, lambda, m, covfn) {
  not_beta <- covfn$not_beta
  not_z <- 1:(4 + covfn$n0)
  beta <- .betaz(pars[-not_beta], id[-not_beta])
  J <- rbind(diag(exp(pars[-not_beta][1:2])), 0)
  J <- cbind(J[, 1], c(0, 0, 1), c(0, 0, 1), J[, 2])
  J <- list(diag(length(not_beta)), J, lapply(sapply(beta[-c(1, 2)], length), diag))
  J <- .blockdiag(.clist(J))
  pars <- c(pars[not_beta], unlist(beta))
  np <- length(pars)
  gH <- .dSigma(pars, XS$X, 0, id, XS$D, covfn, nderiv = 2)
  iAdA <- lapply(gH$d1, function(x) .solve_chol(gH$cholA, x))
  H <- diag(rep(0, np))
  for (i in 1:np) for (j in i:np) {
    temph <- .solve_chol(gH$cholA, gH$d2[[i]][[j]])
    H1 <- (m - 1) * (sum(diag(temph)) - sum(t(iAdA[[i]]) * iAdA[[j]]))
    H2 <- m * sum(diag(.solve_chol(gH$cholA, crossprod(V, 2 * iAdA[[i]] %*% iAdA[[j]] - temph))))
    H[j, i] <- H[i, j] <- .5 * (H1 + H2)
  }
  H <- H + .blockdiag(.clist(list(diag(numeric(4 + covfn$n0)), Map('*', lambda, XS$S))))
  H <- tcrossprod(J %*% H, J)
  list(H = H, J = J, id = id, Vp = MASS::ginv(H))
}

.predict.expansion <- function(object, newx, se.fit = FALSE) {
  if (is.null(newx)) {
    out <- object$fitted
    E <- object$E
    newx <- object$xy
  } else {
    if (!is.null(object$scaling)) {
      newx <- .std_fn(newx, object$scaling)
    }
    E <- .makeE(object$xy, newx)
    out <- newx %*% cbind(object$beta[[1]][1:2], object$beta[[2]][1:2])
    out <- cbind(out, E %*% object$delta)
  }
  if (se.fit) {
    out <- list(fitted = out)
    Vp <- crossprod(object$V$J, object$V$Vp %*% object$V$J)[-object$not_beta, -object$not_beta]
    X <- lapply(object$XS$Uk, function(x) E %*% x)
    X <- .clist(list(newx, newx, X))
    id <- rep(seq_along(X), sapply(X, ncol))
    Vp <- lapply(seq_along(X), function(i) Vp[id == i, id == i])
    ese <- mapply(function(x, y) rowSums(x * t(tcrossprod(y, x))), X, Vp)
    out$se.fit <- sqrt(ese)
  }
  out
}

.plot.expansion <- function(x, start, ptype, brks, pal, onepage, nx, ny, xp, yp, xlab, ylab, ...) {
  if (is.null(xp))
    xp <- pretty(x$xy0[, 1], nx)
  if (is.null(yp))
    yp <- pretty(x$xy0[, 2], ny)
  if (is.null(xlab))
    xlab <- 'x*'
  if (is.null(ylab))
    ylab <- 'y*'
  null.brks <- is.null(brks)
  onebreak <- FALSE
  if (!null.brks & length(brks) == 1) {
    onebreak <- TRUE
    nbreak <- brks
  }
  xyp <- as.matrix(expand.grid(xp, yp))
  xyp2 <- .predict.expansion(x, xyp)
  plots <- lapply(seq_len(ncol(xyp2)), function(i) matrix(xyp2[, i], length(xp)))
  plot_ids <- start:length(plots)
  if (is.null(brks)) {
    brks_lst <- lapply(plots, pretty, n = 10)
  } else {
    if (length(brks) == 1) {
      brks_lst <- lapply(plots, pretty, n = brks)
    } else {
      if (is.vector(brks)) {
        brks_lst <- lapply(seq_along(plot_ids), function(i) brks)
      } else {
        if (is.list(brks)) {
          brks_lst <- brks
        } else {
          stop("'breaks' must be either NULL, an integer, a vector, or a list")
        }
      }
    }
  }
  pal_lst <- lapply(brks_lst, function(x) pal(length(x) - 1))
  if (ptype == 'lattice') {
    grid <- expand.grid(x = xp, y = yp)
    if (!onepage) {
      plot_i <- 0
      for (i in seq_along(plot_ids)) {
        plot_i <- plot_i + 1
        grid$z <- as.vector(plots[[plot_ids[i]]])
        print(lattice::levelplot(z ~ x * y, grid, 
                                 at = brks_lst[[i]], col.regions = pal_lst[[i]]))
      }
    } else { # multiple pages
      lattices <- list()
      for (i in seq_along(plot_ids)) {
        grid$z <- as.vector(plots[[plot_ids[i]]])
        lattices[[i]] <- lattice::levelplot(z ~ x * y, grid, 
                                 at = brks_lst[[i]], col.regions = pal_lst[[i]])
      }
      nxy <- n2mfrow(length(plot_ids))
      do.call(gridExtra::grid.arrange, c(lattices, nrow = nxy[1], ncol = nxy[2]))
    }
  } else {
    for (i in plot_ids) {
      ttl <- substitute(g[i](x), list(i = i))
      image(xp, yp, plots[[i]], xlab = xlab, ylab = ylab, 
            breaks = brks_lst[[i]], col = pal_lst[[i]], ...)
      title(ttl)
    }
  }
}

.search_expansion <- function(pars, id, XS, V, lambda, m, covfn) {
gH <- .d12_expansion(pars, id, XS, V, lambda, m, covfn)
out <- .solve_chol2(gH$H, gH$g)
attr(out, 'gradient') <- gH[[1]]
attr(out, 'Hessian') <- gH[[2]]
out
}

.beta2delta_expansion <- function(beta, UZl) {
sapply(1:2, function(i) as.vector(UZl[[i]] %*% beta[[3]][[i + 1]]))
}
