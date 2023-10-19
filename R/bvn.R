.nllh_cubvn <- function(pars, x, y, ux, uy, freq) {
  pars <- 2 / (1 + exp(-pars)) - 1
  out <- .nllh_bvn_censored_ogram(c(0, pars), x, y, ux, uy, freq)
  if (!is.finite(out)) return(1e20)
  if (out > 1e20) return(1e20)
  out
}

.dcnorm <- function(x, u, mu = 0, sigma = 1, log = FALSE) {
  out <- numeric(length(x))
  above <- x > u
  out[above] <- dnorm(x[above], mu, sigma, log = log)
  out[!above] <- pnorm(u[!above], mu, sigma, log.p = log)
  out
}

.nllh_dcnorm <- function(lsig, x, u) {
  okay <- is.finite(x) & is.finite(u)
  x <- x[okay]
  u <- u[okay]
  out <- -sum(.dcnorm(x, u, 0, exp(lsig), log = TRUE), na.rm = TRUE)
  if (!is.finite(out))
    out <- 1e20
  out
}

.fit_dcnorm <- function(x, u) {
  exp(nlm(.nllh_dcnorm, log(sd(x, na.rm = TRUE)), x = x, u = u)$estimate)
}
