# user-accessible functions
# - deform
# - aniso
# - expand
# - predict.deform
# - plot.deform
# - variogram
# - simulate.deform
# - cencor
# - cencov
#
#' Fitting low-rank nonstationary spatial Gaussian process models through spatial deformation
#'
#' Function \code{deform} fits a 2-dimensional deformation model, where typically
#' x and y coordinates in geographic (G-) space will be provided and then deformed
#' to give new coordinates in deformed (D-) space in which isotropy of a Gaussian 
#' process is optimally achieved.
#'
#' @param x a 2-column matrix comprising x and y coordinates column-wise, respectively, or a list; see Details for the latter
#' @param z a variance-covariance matrix
#' @param n an integer number of data
#' @param k an integer vector of ranks
#' @param lambda specified lambda values; see Details
#' @param lambda0 initial lambda values
#' @param correlation a logical defining whether \code{z} should be assumed to be a correlation matrix; defaults to \code{FALSE}
#' @param cosine a logical defining whether the powered exponential covariance function should be multiplied by the cosine of scaled distances, i.e. giving a damped oscillation; defaults to \code{FALSE}
#' @param bijective a logical for whether a bijective deformation should be imposed; defaults to FALSE
#' @param bijective.args a list specifying quantities to ensure bijectivity, if bijective == TRUE; see Details
#' @param trace an integer specifying the amount to report on optimisation (0, default, is nothing; 1 gives a bit)
#' @param standardise a character string that governs whether dimensions are scaled by a common (\code{"together"}) or dimension-specific factor; defaults to \code{"together"}
#'
#' @details
#'
#' If \code{x} is a list, then it wants elements \code{"x"}, \code{"z"} and \code{"n"} as described above.
#' 
#' Values of \code{lambda} multiply the penalties placed on the wiggliness of the 
#' smooths that form the deformations. Larger values make things less wiggly. Values 
#' of \code{lambda0} specify initial values for \code{lambda}, which are still optimised.
#'
#' \code{bijective.args()} is a 4-element list: \code{"mult"} is a penalty placed on
#' the numerical approximation to identifying non-bijectivity, where larger values
#' impose bijectivity more strictly; \code{"scl"} is a scaling placed on the grid
#' used to numerically identify non-bijectivity, where smaller values will typically 
#' impose bijectivity more strictly; \code{"nx"} and \code{"ny"} specify the x and y
#' dimensions of the grid used to numerically identify bijectivity. Defaults are 
#' \code{mult = 1e3}, \code{scl = 1}, \code{nx = 40} and \code{ny = 40}. It is advisable
#' to use \code{"mult"} and not \code{"scl"} to control bijectivity, in the first instance.
#' 
#' @references
#'
#' Sampson, P. D. and Guttorp, P. (1992) Nonparametric Estimation of Nonstationary Spatial
#' Covariance Structure, Journal of the American Statistical Association, 87:417, 108-119,
#' \doi{10.1080/01621459.1992.10475181}
#' 
#' Wood, S.N. (2003), Thin plate regression splines. Journal of the Royal Statistical 
#' Society: Series B (Statistical Methodology), 65: 95-114. 
#' \doi{10.1111/1467-9868.00374}
#' 
#' @examples
#'
#' \donttest{
#' 
#' data(solar)
#' deform(solar$x, solar$z, solar$n) 
#' # equivalent to deform(solar)
#' 
#' # bijective deformation 
#' deform(solar, bijective = TRUE)
#' 
#' # deformation with specified rank
#' deform(solar, k = c(10, 8))
#' 
#' }
#'
#' @return An object of class \code{deform}  and then of class \code{deformation}
#'
#' @export
#'
deform <- function(x, z, n, k = c(10, 10), lambda = c(-1, -1), 
lambda0 = rep(exp(3), length(k)), correlation = FALSE, cosine = FALSE, 
bijective = FALSE, bijective.args = NULL, trace = 0, standardise = 'together') {
if (is.list(x)) {
  z <- x$z
  n <- x$n
  x <- x$x
}
if (max(k) > nrow(x))
  stop("Can't have a higher low-rank representation than the number of data points!")
x0 <- x
x <- .std_fn(x, NULL, standardise == 'together')
M <- 3
E <- .makeE(x)
eE <- eigen(E, symmetric = TRUE)
T <- cbind(1, x)
XS0 <- .makeXS(eE, T, k, M)
id <- rep(seq_len(length(k) + 1), c(M, k - M))
.cov_fns <- .make_cov_fns(correlation, cosine)
inits <- .inits(x[, 1], x[, 2], z, n, .cov_fns)
inits <- c(inits, numeric(sum(k - M) + 1))
XS0$D <- lapply(XS0$X, function(z) lapply(seq_len(ncol(z)), function(i) .vecdiff(z[, i])))
ctrl <- list(steptol = 1e-12, itlim = 1e2, fntol = 1e-8, gradtol = 1e-2, stepmax = 3)
fit_lambda <- lambda < 0
if (bijective)
  attr(bijective, 'data') <- .bijective_grid(x, XS0, inits, bijective.args)
if (all(fit_lambda)) {
  rho0 <- log(lambda0)
  attr(rho0, 'beta') <- inits
  f1 <- .BFGS(rho0, .reml0, .reml1, id = id, XS = XS0, V = z, m = n, 
              covfn = .cov_fns, biject = bijective, control = ctrl, trace = trace)
  bb <- attr(f1$objective, 'beta')
  lambda <- exp(as.vector(f1$par))
} else {
  if (all(!fit_lambda)) {
    f1 <- .newton_step_inner(inits, .d0_deform, .search_deform, id = id, XS = XS0, 
                             V = z, lambda = lambda, m = n, covfn = .cov_fns, 
                              biject = bijective)
    bb <- f1$par
  } else {
    rho0 <- log(lambda0[fit_lambda])
    attr(rho0, 'beta') <- inits
    attr(lambda, 'fit') <- fit_lambda
    f1 <- .BFGS(rho0, .reml0, .reml1, id = id, XS = XS0, V = z, m = n, 
                lambda = lambda, covfn = .cov_fns, biject = bijective, 
                control = ctrl, trace = trace)
    lambda[fit_lambda] <- exp(as.vector(f1$par))
    bb <- attr(f1$objective, 'beta')
  }
}
beta <- .betaxy(bb[-.cov_fns$not_beta], id)
delta <- .beta2delta(beta, XS0$UZ)
xy_hat <- .new_xy(bb[-.cov_fns$not_beta], id, XS0$X, atts = FALSE)
cov_pars <- bb[.cov_fns$not_beta]
cov_pars0 <- c(p0 = 1e3, p1 = 0, p2 = NA, p3 = NA)
cov_pars <- replace(cov_pars0, .cov_fns$nms[.cov_fns$not_beta], bb[.cov_fns$not_beta])
cov_pars <- c(exp(cov_pars[1]), exp(cov_pars[2]), .9996 / (1 + exp(-cov_pars[3])), 2 / (1 + exp(-cov_pars[4])))
out <- list(beta = beta, XS = XS0, delta = delta, xy = x, fitted = xy_hat, z = z, 
            cov_pars = cov_pars, lambda = lambda, xy0 = x0, E = E,
            covfn = .cov_fns, cosine = cosine, correlation = correlation,
            not_beta = .cov_fns$not_beta, b0 = bb[.cov_fns$not_beta], 
            scaling = attr(x, 'scaling'))
out$V <- .vcov_deformation(bb, id = id, XS = XS0, V = z, m = n, lambda = lambda, 
                           covfn = .cov_fns, biject = bijective)
class(out) <- c('deform', 'deformation')
out
}

#' Fitting anisotropic spatial Gaussian process models
#'
#' Function \code{aniso} fits a conventional 2-dimensional anisotropic Gaussian 
#' process, i.e. just with scalings in the x and y coordinates.
#'
#' @param x a 2-column matrix comprising x and y coordinates column-wise, respectively, or a list; see Details for the latter
#' @param z a variance-covariance matrix
#' @param n an integer number of data
#' @param correlation a logical defining whether \code{z} should be assumed to be a correlation matrix; defaults to \code{FALSE}
#' @param cosine a logical defining whether the powered exponential covariance function should be multiplied by the cosine of scaled distances, i.e. giving a damped oscillation; defaults to \code{FALSE}
#' @param standardise a character string that governs whether dimensions are scaled by a common (\code{"together"}) or dimension-specific factor; defaults to \code{"together"}
#'
#' @details
#'
#' If \code{x} is a list, then it wants elements \code{"x"}, \code{"z"} and \code{"n"} as described above.
#'
#' @references
#'
#' Sampson, P. D. and Guttorp, P. (1992) Nonparametric Estimation of Nonstationary Spatial
#' Covariance Structure, Journal of the American Statistical Association, 87:417, 108-119,
#' \doi{10.1080/01621459.1992.10475181}'
#' 
#' @examples
#'
#' data(solar)
#' aniso(solar$x, solar$z, solar$n) 
#' # equivalent to aniso(solar)
#' 
#' @return An object of class \code{deform} and then of class \code{anisotropic}
#'
#' @export
#'
aniso <- function(x, z, n, correlation = FALSE, cosine = FALSE, standardise = 'together') {
if (is.list(x)) {
  z <- x$z
  n <- x$n
  x <- x$x
}
x <- .std_fn(x, NULL, standardise == 'together')
.cov_fns <- .make_cov_fns(correlation, cosine)
bb <- .inits(x[, 1], x[, 2], z, n, .cov_fns, hessian = TRUE)
V <- MASS::ginv(attr(bb, 'hessian'))
not_beta <- .cov_fns$not_beta
bb[-not_beta] <- exp(bb[-not_beta])
xy_hat <- t(t(x) * bb[-not_beta])
cov_pars <- bb[.cov_fns$not_beta]
cov_pars0 <- c(p0 = 1e3, p1 = 0, p2 = NA, p3 = NA)
cov_pars <- replace(cov_pars0, .cov_fns$nms[.cov_fns$not_beta], bb[.cov_fns$not_beta])
cov_pars <- c(exp(cov_pars[1]), exp(cov_pars[2]), 1 / (1 + exp(-cov_pars[3])), 2 / (1 + exp(-cov_pars[4])))
out <- list(xy = x, fitted = xy_hat, z = z, cov_pars = cov_pars, b = bb[-not_beta],
            covfn = .cov_fns, cosine = cosine, correlation = correlation,
            not_beta = not_beta, b0 = bb[not_beta], scaling = attr(x, 'scaling'))
out$V <- list(Vp = V)
class(out) <- c('deform', 'anisotropic')
out
}

#' Fitting low-rank nonstationary spatial Gaussian process models through dimension expansion
#'
#' Function \code{exapnd} fits a multi-dimensional dimension expansion model, where typically
#' x and y coordinates in geographic (G-) space will be provided and then scaled and
#' combined with new latent dimensions (that a functions of x and y) to give new coordinates 
#' in deformed (D-) space in which isotropy of a Gaussian process is optimally achieved.
#'
#' @param x a 2-column matrix comprising x and y coordinates column-wise, respectively, or a list; see Details for the latter
#' @param z a variance-covariance matrix
#' @param n an integer number of data
#' @param k an integer vector of ranks
#' @param lambda specified lambda values
#' @param lambda0 initial lambda values
#' @param correlation a logical defining whether \code{z} should be assumed to be a correlation matrix; defaults to \code{FALSE}
#' @param cosine a logical defining whether the powered exponential covariance function should be multiplied by the cosine of scaled distances, i.e. giving a damped oscillation; defaults to \code{FALSE}
#' @param trace an integer specifying the amount to report on optimisation (0, default, is nothing; 1 gives a bit)
#' @param z0 a scalar giving initial values (which alternate \code{z0, -z0, z0, ...} for latent dimensions
#' @param standardise a character string that governs whether dimensions are scaled by a common (\code{"together"}) or dimension-specific factor; defaults to \code{"together"}
#'
#' @details
#'
#' If \code{x} is a list, then it wants elements \code{"x"}, \code{"z"} and \code{"n"} as described above.
#'
#' @references
#'
#' Bornn, L., Shaddick, G., & Zidek, J. V. (2012). Modeling nonstationary processes through
#' dimension expansion. Journal of the American Statistical Association, 107(497), 281-289.
#' \doi{10.1080/01621459.2011.646919}.
#' 
#' @examples
#' 
#' # one-dimensional expansion
#' data(solar)
#' expand(solar$x, solar$z, solar$n)
#' # equivalent to expand(solar)
#' 
#' \donttest{
#' 
#' # two-dimensional expansion with rank-8 and rank-5 dimensions
#' expand(solar$x, solar$z, solar$n, c(8, 5))
#' 
#' }
#' 
#' @return An object of class \code{deform} and then of class \code{expansion}
#'
#' @export
#'
expand <- function(x, z, n, k = 10, lambda = rep(-1, length(k)), 
                   lambda0 = rep(exp(3), length(k)), correlation = FALSE, 
                   cosine = FALSE, trace = 0, z0 = NULL, standardise = 'together') {
  if (is.list(x)) {
    z <- x$z
    n <- x$n
    x <- x$x
  }
  if (max(k) > nrow(x))
    stop("Can't have a higher low-rank representation than the number of data points!")
  x0 <- x
  x <- .std_fn(x, NULL, standardise == 'together')
  E <- .makeE(x)
  eE <- eigen(E, symmetric = TRUE)
  XS0 <- .makeXS_expansion(eE, k, x)
  .cov_fns <- .make_cov_fns(correlation, cosine)
  id <- c(1:.cov_fns$n0, rep(.cov_fns$n0 + 1:(length(k) + 2), c(2, 2, k)))
  inits <- .inits(x[, 1], x[, 2], z, n, .cov_fns)
  if (is.null(z0))
    z0 <- .05 * min(apply(x, 2, function(x) diff(range(x))))
  inits <- c(inits, 0, rep(z0, sum(k)) * rep(c(-1, 1), ceiling(sum(k) / 2))[1:sum(k)])
  XS0$D <- lapply(XS0$X, function(z) lapply(seq_len(ncol(z)), function(i) .vecdiff(z[, i])))
  ctrl <- list(steptol = 1e-12, itlim = 1e2, fntol = 1e-8, gradtol = 1e-2, stepmax = 3)
  fit_lambda <- lambda < 0
  if (all(fit_lambda)) {
    rho0 <- log(lambda0)
    attr(rho0, 'beta') <- inits
    f1 <- .BFGS(rho0, .reml0_expansion, .reml1_expansion, id = id, XS = XS0, V = z, 
                m = n, covfn = .cov_fns, control = ctrl, trace = trace)
    bb <- attr(f1$objective, 'beta')
    lambda <- exp(as.vector(f1$par))
  } else {
    if (all(!fit_lambda)) {
      f1 <- .newton_step_inner(inits, .d0_expansion, .search_expansion, id = id, 
                               XS = XS0, V = z, lambda = lambda, m = n, covfn = .cov_fns)
      bb <- f1$par
    } else {
      rho0 <- log(lambda0[fit_lambda])
      attr(rho0, 'beta') <- inits
      attr(lambda, 'fit') <- fit_lambda
      f1 <- .BFGS(rho0, .reml0_expansion, .reml1_expansion, id = id, XS = XS0, V = z, 
                  m = n, covfn = .cov_fns, lambda = lambda, control = ctrl, trace = trace)
      lambda[fit_lambda] <- exp(as.vector(f1$par))
      bb <- attr(f1$objective, 'beta')
    }
  }
  beta <- .betaz(bb[-.cov_fns$not_beta], id[-.cov_fns$not_beta])
  xyz_hat <- mapply('%*%', XS0$X, beta)
  delta <- mapply('%*%', XS0$U, beta[-c(1, 2)])
  cov_pars <- bb[.cov_fns$not_beta]
  cov_pars0 <- c(p0 = 1e3, p1 = 0, p2 = NA, p3 = NA)
  cov_pars <- replace(cov_pars0, .cov_fns$nms[.cov_fns$not_beta], bb[.cov_fns$not_beta])
  cov_pars <- c(exp(cov_pars[1]), exp(cov_pars[2]), 1 / (1 + exp(-cov_pars[3])), 2 / (1 + exp(-cov_pars[4])))
  out <- list(beta = beta, XS = XS0, xy = x, fitted = xyz_hat, delta = delta, z = z, 
              cov_pars = cov_pars, lambda = lambda, xy0 = x0, E = E,
              covfn = .cov_fns, cosine = cosine, correlation = correlation,
              not_beta = .cov_fns$not_beta, b0 = bb[.cov_fns$not_beta], 
              scaling = attr(x, 'scaling'))
  out$V <- .vcov_expansion(bb, id = id, XS = XS0, V = z, lambda = lambda, m = n,
                           covfn = .cov_fns)
  #out$scaling <- scl
  class(out) <- c('deform', 'expansion')
  out
}

#' Predict from a fitted \code{deform} object
#'
#' @param object a fitted \code{deform} object
#' @param newdata a 2-column matrix of x and y coordinates
#' @param ... currently just a placeholder
#'
#' @return A 2-column matrix of predicted x and y points for deformations and
#' a (2 + q)-column matrix for q-dimensional expansions.
#'
#' @examples
#'
#' \donttest{
#'
#' # fit a deformation model
#' data(solar)
#' m0 <- deform(solar$x, solar$z, solar$n)
#' 
#' # predict D-space points for original locations
#' predict(m0)
#' 
#' }
#' 
#' # predictions for one-dimensional expansion model with specified locations 
#' # and standard error estimates
#' data(solar)
#' m1 <- expand(solar$x, solar$z, solar$n)
#' xvals <- seq(-123.3, -122.2, by = .1)
#' yvals <- seq(49, 49.4, by = .1)
#' xyvals <- expand.grid(xvals, yvals)
#' predict(m1, xyvals, se.fit = TRUE)
#'
#' @export
#'
predict.deform <- function(object, newdata = NULL, ...) {
cls <- class(object)
if (cls[2] == 'deformation') {
  out <- .predict.deformation(object, newdata, ...)
} else {
  if (cls[2] == 'anisotropic') {
    out <- .predict.aniso(object, newdata, ...)
  } else {
    if (cls[2] == 'expansion') {
      out <- .predict.expansion(object, newdata, ...)
    }
  }
}
out
}

#' Plot a fitted \code{deform} object
#'
#' @param x a fitted \code{deform} object
#' @param start an integer giving the starting dimension of plots of dimension expansion models; defaults to 1
#' @param graphics a character string that is either \code{"graphics"} or \code{"lattice"} and states the graphics package to use for plots; defaults to \code{"graphics"}
#' @param breaks an integer, vector or list; see Details
#' @param pal a function specifying the colour palette to use for plots; defaults to \code{hcl.colors(..., 'YlOrRd', rev = TRUE)}
#' @param onepage a logical specifying whether all plots should be put on one page; defaults to \code{FALSE}, which makes use of the current graphics state
#' @param nx number of x points to use for plotting grid
#' @param ny number of y points to use for plotting grid
#' @param xp x points to use for plotting grid
#' @param yp y points to use for plotting grid
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ... extra arguments to pass to \code{plot()}
#' 
#' @details
#' 
#' If \code{breaks} is an integer then it specifies the number of breaks to use for colour scales; if it's a vector, then it's the breaks themselves; and if it's a list then it's different breaks for each dimension.
#'
#' @return Plots representing all one- or two-dimensional smooths
#'
#' @examples
#'
#' \donttest{
#' 
#' # deformations
#' data(solar)
#' m0 <- deform(solar$x, solar$z, solar$n)
#' 
#' # plot representation of deformation
#' plot(m0)
#' 
#' # as above with specified x and y grid
#' xvals <- seq(-123.3, -122.25, by = .05)
#' yvals <- seq(49, 49.4, by = .05)
#' plot(m0, xp = xvals, yp = yvals)
#' 
#' }
#' 
#' # one-dimensional expansion
#' data(solar)
#' m1 <- expand(solar$x, solar$z, solar$n)
#' 
#' # plot its three dimensions
#' op <- par(mfrow = c(1, 3))
#' plot(m1)
#' par(op)
#' 
#' # or plot using lattice::levelplot
#' plot(m1, graphics = 'lattice')
#' # or as above, but on one page
#' plot(m1, graphics = 'lattice', onepage = TRUE)
#' 
#' \donttest{
#' 
#' # two-dimensional expansion
#' m2 <- expand(solar$x, solar$z, solar$n, c(8, 5)) 
#' # plot of its third and fourth dimensions for given x and y values
#' op <- par(mfrow = c(1, 2))
#' plot(m2, start = 3, xp = xvals, yp = yvals)
#' par(op)
#' 
#' # using lattice::levelplot with common breaks across dimensions with
#' # a palette that gives latent dimensions in white where near zero
#' plot(m2, onepage = TRUE, graphics = 'lattice', breaks = seq(-0.35, 0.35, by = 0.1), 
#'      pal = function(n) hcl.colors(n, 'Blue-Red 3'))
#'      
#' }
#' 
#' @export
#'
plot.deform <- function(x, start = 1, graphics = 'base', breaks = NULL, 
                        pal = function(n) hcl.colors(n, 'YlOrRd', rev = TRUE), 
                        onepage = FALSE, nx = 10, ny = 10, xp = NULL, yp = NULL, 
                        xlab = NULL, ylab = NULL, ...) {
cls <- class(x)
if (cls[2] == 'deformation') .plot.deformation(x, nx, ny, xp, yp, xlab, ylab, ...)
if (cls[2] == 'anisotropic') .plot.aniso(x, nx, ny, xp, yp, xlab, ylab, ...)
if (cls[2] == 'expansion') .plot.expansion(x, start, graphics, 
                                           breaks, pal, onepage, nx, ny, xp, yp, xlab, ylab, ...)
}

#' Plot the variogram for a fitted \code{deform} object
#'
#' @param object a fitted \code{deform} object
#' @param bins an integer specifying the number of bins for plotting
#' @param bin.function a character specifying a function to use to calculate bins; defaults to \code{pretty()}
#' @param trim a scalar in [0, 0.5], which is passed to \code{mean()} when calculating binned variogram estimates; defaults to 0
#' @param ... extra arguments to pass to \code{plot()}
#' 
#'
#' @return Plot of variogram
#'
#' @examples
#' 
#' \donttest{
#' 
#' # deformations
#' data(solar)
#' m0 <- deform(solar$x, solar$z, solar$n)
#'
#' # empirical versus model-based variogram estimates against distance,
#' # where distance is based on D-space
#' variogram(m0)
#' # which is the default with approximately 20 bins, i.e. variogram(m0, bins = 20)
#' 
#' }
#' 
#' # variogram for one-dimensional expansion without binning
#' data(solar)
#' m1 <- expand(solar$x, solar$z, solar$n)
#' variogram(m1, bins = 0)
#'
#' @export
#'
variogram <- function(object, bins = 20, bin.function = 'pretty', trim = 0, ...) {
xy <- object$fitted
pars <- object$cov_pars
sv <- .semivariog(xy, object$z, bins = bins, bin.function, trim)
plot(sv[, 1], sv[, 2], xlab = 'Distance', ylab = 'Semivariance', ylim = c(0, max(sv[, 2], na.rm = TRUE)), ...)
xp <- pretty(par('usr')[1:2], 100)
if (object$cosine) {
  lines(xp, pars[2] * (1 - pars[3] * exp(-xp^pars[4]) * cos(xp / pars[1])))
} else {
  lines(xp, pars[2] * (1 - pars[3] * exp(-xp^pars[4])))
}    
}

#' Simulate from a fitted \code{deform} object
#'
#' @param object a fitted \code{deform} object
#' @param nsim an integer giving the number of simulations
#' @param seed an integer giving the seed for simulations
#' @param newdata a 2-column matrix of x and y coordinates
#' @param ... extra arguments to pass to \code{predict.deform()}
#'
#' @return Plots representing all one- or two-dimensional smooths
#'
#' @examples
#' 
#' \donttest{
#'
#' # deformations
#' data(solar)
#' m0 <- deform(solar$x, solar$z, solar$n) 
#' # Gaussian process simulations based on fitted deformation model
#' simulate(m0)
#' 
#' }
#' 
#' # one-dimensional expansion model with five simulations and specified locations
#' data(solar)
#' m1 <- expand(solar$x, solar$z, solar$n)
#' xvals <- seq(-123.3, -122.25, by = .05)
#' yvals <- seq(49, 49.4, by = .05)
#' xyvals <- expand.grid(xvals, yvals)
#' simulate(m1, 5, newdata = xyvals)
#'
#' @export
#'
simulate.deform <- function(object, nsim = 1, seed = NULL, newdata = NULL, ...) {
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1) # initialize the RNG if necessary
  if(is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  xy <- predict.deform(object, newdata, ...)
  V <- object$covfn$pd0(object$b0, xy, nrow(xy))
  .pivchol_rmvn(nsim, numeric(nrow(V)), V)
}

#' Correlation and covariance matrices from censored data
#'
#' @param x a numeric matrix
#' @param u a numeric matrix giving corresponding points of left-censoring
#'
#' @return a matrix
#' 
#' @details For \code{cencov()} a covariance matrix is returned and for 
#' \code{cencor()} a correlation matrix is returned. Note that \code{cencov()} 
#' calls \code{cencor()}. Estimates are based on assuming values are from a 
#' multivariate Gaussian distribution.
#' 
#' @seealso \link{cov} and \link{cor} for uncensored estimates.
#'
#' @examples
#'
#' # generate some correlated data
#' n <- 1e2
#' x <- rnorm(n)
#' y <- 0.25 * x + sqrt(0.75) * rnorm(n)
#' xy <- cbind(x, y)
#' # threshold of zero for left-censoring
#' u <- matrix(0, n, 2) 
#' # left-censored correlation matrix
#' cencor(xy, u) # could check with cor(xy)
#' # left-censored covariance matrix
#' cencov(xy, u)
#'
#' @export
#'
cencov <- function(x, u) {
  u <- as.matrix(u)
  if (ncol(u) == 1 & ncol(x) != 1) {
    u <- matrix(u, nrow(x), ncol(x))
  }
  out <- sapply(seq_len(ncol(x)), function(i) .fit_dcnorm(x[, i], u[, i]))
  x <- x / out
  u <- u / out
  out <- tcrossprod(out)
  out * cencor(x, u)
}

#'
#' @rdname cencov
#' 
#' @export
#' 
cencor <- function(x, u) {
  u <- as.matrix(u)
  if (ncol(u) == 1 & ncol(x) != 1) {
    u <- matrix(u, nrow(x), ncol(x))
  }
  out <- diag(1, ncol(x))
  for (i in 1:ncol(x)) for (j in 1:i) {
    if (j < i) {
      temp <- nlm(function(z) .nllh_cubvn(z, x[, i], x[, j], u[, i], u[, j], 1), .5)$estimate
      out[i, j] <- 2 / (1 + exp(-temp)) - 1
      out[j, i] <- out[i, j]
    }
  }
  out
}
