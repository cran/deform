# .newton_step_inner
# .cholPDchol
# .logdetH
# .solve_chol
# .itreport
# .BFGS

.newton_step_inner <- function(pars, fn, sfn, ..., steptol=1e-12, itlim=1e2, fntol=1e-8,
gradtol=1e-4, trace=0, stepmax=1e8) {

pars0 <- pars

it <- 1
okay <- TRUE
f0 <- fn(pars, ...)

while (okay) {
if (it > 1) g0 <- g
if (exists("step1")) {
step0 <- step1
g <- attr(step0, "gradient")
} else {
step0 <- sfn(pars, ...)
g <- attr(step0, "gradient")
}
if (mean(abs(g)) < gradtol) {
report <- c("gradient tolerance reached")
break
}

step0 <- sign(step0) * pmin(abs(step0), stepmax)
alpha <- 1
report <- NULL
ls <- TRUE
while(ls & is.null(report)) {
step <- alpha * step0
stepokay <- all(abs(step) > steptol)
if (!stepokay) {
report <- c("step tolerance reached")
} else {
theta1 <- pars - step
f1 <- fn(theta1, ...)
d <- f1 - f0
if (d < 0) {
step1 <- sfn(theta1, ...)
if (any(!is.finite(attr(step1, "gradient")))) d <- 1
}
if (d < 0) {
f0 <- f1
pars <- theta1
ls <- FALSE
} else {
if (d < fntol) {
report <- c("function tolerance reached")
}
alpha <- .5 * alpha
}
}
}
if (!is.null(report)) break
it <- it + 1
if (it == itlim) {
report <- c("iteration limit reached")
okay <- FALSE
}
}
g <- attr(step0, "gradient")
H <- attr(step0, "Hessian")
iH <- attr(step0, "invHessian")
list(par=as.vector(pars), objective=f0, gradient=g, hessian=H, report=report, convergence=0, iterations=it, beta=attr(pars, "beta"), invHessian=iH)
}

.precondition <- function(H) {
d <- 1 / sqrt(diag(abs(H)))
D <- diag(d)
H <- D %*% H %*% D
attr(H, 'd') <- d
H
}

.cholPDchol <- function(A) {
d0 <- diag(A)
eps <- 1e-16
test <- try(chol(A), silent=TRUE)
while(inherits(test, "try-error")) {
diag(A) <- d0 + eps
test <- try(chol(A), silent=TRUE)
eps <- 2 * eps
}
test
}

.cholPDchol <- function(A) {
d0 <- diag(A)
eps <- 1e-16
test <- try(chol(A), silent=TRUE)
while(inherits(test, "try-error")) {
diag(A) <- d0 + eps
test <- try(chol(A), silent=TRUE)
eps <- 2 * eps
}
test
}

.logdetH <- function(H) {
D <- diag(1 / sqrt(abs(diag(H))))
cholH <- .cholPDchol(D %*% H %*% D)
ldetH <- 2 * sum(log(diag(cholH) / diag(D)))
out <- list(ldetH=ldetH)
return(out)
}

.solve_chol <- function(L, x) {
piv <- ipiv <- attr(L, "pivot")
if (is.null(piv)) {
  out <- backsolve(L, backsolve(L, x, transpose = TRUE))
} else {
  ipiv[piv] <- seq_along(piv)
  out <- backsolve(L, backsolve(L, x[piv, , drop=FALSE], upper.tri=TRUE, transpose=TRUE))[ipiv, , drop=FALSE]
}
out
}

.solve_pchol_d <- function(L, d, x) {
D <- diag(d)
x <- d * x
piv <- ipiv <- attr(L, "pivot")
ipiv[piv] <- seq_along(piv)
out <- backsolve(L, backsolve(L, x[piv, , drop=FALSE], upper.tri=TRUE, transpose=TRUE))[ipiv, , drop = FALSE]
D %*% out
}

.itreport <- function(f, g, it) {
    report <- paste("\n Outer iteration ", it, ":", sep="")
    rep1 <- paste("  Outer max(|grad|):", signif(max(abs(g)), 3))
    rep2 <- paste("  Inner max(|grad|): ", signif(max(abs(attr(f, "gradient"))), 3), ".", sep="")
    report <- c(report, paste(rep1, rep2, sep="; "))
    cat(paste(report, collapse="\n"))
}

.BFGS <- function(pars, fn, gfn, ..., control, trace=0) {

steptol <- control$steptol
itlim <- control$itlim
fntol <- control$fntol
gradtol <- control$gradtol
stepmax <- control$stepmax

it <- 1
okay <- TRUE
f0 <- fn(pars, ...)
g1 <- NULL
I <- iH <- H <- diag(length(pars))

while (okay) {
if (it > 1) g0 <- g
if (!is.null(g1)) {
g <- g1
} else {
attr(pars, "beta") <- attr(f0, "beta")
g <- gfn(pars, ...)
}
if (trace) .itreport(f0, g, it - 1)
if (mean(abs(g)) < gradtol) {
report <- c("gradient tolerance reached")
break
}
step0 <- crossprod(iH, g)
step0 <- sign(step0) * pmin(abs(step0), stepmax)
alpha <- 1
report <- NULL
ls <- TRUE
while(ls & is.null(report)) {
step <- alpha * step0
stepokay <- all(abs(step) > steptol)
if (!stepokay) {
report <- c("step tolerance reached")
} else {
theta1 <- pars - step
f1 <- fn(theta1, ...)
d <- f1 - f0
if (d < 0) {
attr(theta1, "beta") <- attr(f1, "beta")
g1 <- gfn(theta1, ...)
if (any(!is.finite(g1))) d <- 1
yk <- g1 - g
denom <- sum(- yk * step)
t1 <- I - tcrossprod(- step, yk) / denom
t2 <- I - tcrossprod(yk, - step) / denom
t3 <- tcrossprod(- step) / denom
iH <- t1 %*% iH %*% t2 + t3
if (any(!is.finite(iH))) d <- 1
}
if (d < 0) {
f0 <- f1
pars <- theta1
ls <- FALSE
} else {
if (d < fntol) {
report <- c("function tolerance reached")
}
alpha <- .5 * alpha
}
}
}
if (!is.null(report)) break
it <- it + 1
if (it == itlim) {
report <- c("iteration limit reached")
okay <- FALSE
}
}
if (trace) cat(paste("\n ", it, "iterations:", report, "\n"))
out <- list(par=as.vector(pars), objective=f0)
out$gradient <- g
out$convergence <- 0
out$report <- report
out$iterations <- it
if (!is.null(attr(pars, "beta"))) out$beta <- attr(pars, "beta")
out
}

.pivchol_rmvn <- function(n, mu, Sig) {
  R <- suppressWarnings(chol(Sig, pivot = TRUE))
  piv <- order(attr(R, "pivot"))  ## reverse pivoting index
  r <- attr(R, "rank")  ## numerical rank
  V <- R[1:r, piv]
  Y <- crossprod(V, matrix(rnorm(n * r), r))
  Y + as.vector(mu)
}