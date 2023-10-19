.bijective_grid <- function(x, XS, inits, lst = NULL) {

lst0 <- list(nx = 40, ny = 40, scl = 1, mult = 1e3)

if (is.null(args)) {
  lst <- lst0
} else {
  lst <- replace(lst0, names(lst), lst)
}

nx <- lst$nx
ny <- lst$ny

xx <- seq(min(x[, 1]) - 1e-6, max(x[, 1]) + 1e-6, l = nx + 1)
yy <- seq(min(x[, 2]) - 1e-6, max(x[, 2]) + 1e-6, l = ny + 1)

newx <- as.matrix(expand.grid(yy, xx))[, 2:1]

hx <- xx[2] - xx[1]
hy <- yy[2] - yy[1]
hxy <- hx * hy
scl <- .5 * hxy * lst$scl * prod(exp(tail(inits[5:6], 2)))

E <- .makeE(x, newx)
X <- lapply(XS$UZ, function(x) cbind(newx, E %*% x))
T <- cbind(1, newx)

xm <- matrix(newx[, 1], nx + 1)
ym <- matrix(newx[, 2], nx + 1)
idm <- matrix(seq_len(nrow(newx)), nx + 1, byrow = TRUE)

ij <- ij0 <- list()
it <- 1
for (i in seq_len(nx)) {
  for (j in seq_len(ny)) {
    temp <- idm[cbind(i + c(0, 0, 1, 1), j + c(0, 1, 1, 0))]
    ij0[[it]] <- temp
    ij[[it]] <- temp[c(1:3, 1, 3:4)]
    it <- it + 1
  }
}

nvec <- rep(seq_along(X), sapply(X, ncol))

out <- list(x = newx, hxy = hxy, E = E, T = T, ij = unlist(ij), mult = lst$mult, X = X, scl = scl, nvec = nvec)

}

.area2pen_d0 <- function(x, scl = .5) {
  x <- x / scl
  ifelse(x < 0, x^2, 0)
}

.area2pen_d1 <- function(x, scl = .5) {
  x <- x / scl
  ifelse(x < 0, 2 * x, 0) / scl
}

.area2pen_d2 <- function(x, scl = .5) {
  x <- x / scl
  ifelse(x < 0, 2, 0) / scl^2
}

.pen_deriv <- function(pars, X1, X2, nvec, id, deriv = 0, scl) {
beta <- split(pars, nvec)
x <- X1 %*% beta[[1]]
y <- X2 %*% beta[[2]]
x123 <- matrix(x[id], 3)
y123 <- matrix(y[id], 3)
d0 <- x123[2, ] * y123[1, ] + x123[3, ] * y123[2, ] + x123[1, ] * y123[3, ]
d0 <- d0 - x123[1, ] * y123[2, ] - x123[2, ] * y123[3, ] - x123[3, ] * y123[1, ]
d0 <- .5 * d0
out <- list(d0 = sum(.area2pen_d0(d0, scl)))
if (1 %in% deriv | 2 %in% deriv) {
  idm <- matrix(id, 3)
  d11 <- X1[idm[1, ], ] * (y123[3, ] - y123[2, ])
  d11 <- d11 + X1[idm[2, ], ] * (y123[1, ] - y123[3, ])
  d11 <- d11 + X1[idm[3, ], ] * (y123[2, ] - y123[1, ])
  d12 <- X2[idm[1, ], ] * (x123[2, ] - x123[3, ])
  d12 <- d12 + X2[idm[2, ], ] * (x123[3, ] - x123[1, ])
  d12 <- d12 + X2[idm[3, ], ] * (x123[1, ] - x123[2, ])
  d1 <- .5 * cbind(d11, d12)
  if (1 %in% deriv)
    out$g <- colSums(.area2pen_d1(d0, scl) * d1)
  if (2 %in% deriv) {
    d212 <- crossprod(X1[idm[1, ], ], .area2pen_d1(d0) * (X2[idm[3, ], ] - X2[idm[2, ], ]))
    d212 <- d212 + crossprod(X1[idm[2, ], ], .area2pen_d1(d0) * (X2[idm[1, ], ] - X2[idm[3, ], ]))
    d212 <- d212 + crossprod(X1[idm[3, ], ], .area2pen_d1(d0) * (X2[idm[2, ], ] - X2[idm[1, ], ]))
    d212 <- .5 * d212
    d211 <- diag(numeric(ncol(X1)))
    d222 <- diag(numeric(ncol(X2)))
    d2 <- rbind(cbind(d211, d212), cbind(t(d212), d222))
    out$H <- d2 + crossprod(d1, .area2pen_d2(d0, scl) * d1)
  }
}
out
}

.d0_biject <- function(pars, id, XS, biject) {
bij_data <- attr(biject, 'data')
out <- .pen_deriv(pars, bij_data$X[[1]], bij_data$X[[2]], bij_data$nvec, bij_data$ij, 0, bij_data$scl)$d0
bij_data$mult * out
}

.d12_biject <- function(pars, id, XS, biject) {
bij_data <- attr(biject, 'data')
gH <- .pen_deriv(pars, bij_data$X[[1]], bij_data$X[[2]], bij_data$nvec, bij_data$ij, 1:2, bij_data$scl)
gH$g <- bij_data$mult * gH$g
gH$H <- bij_data$mult * gH$H
gH
}
