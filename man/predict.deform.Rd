% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3.R
\name{predict.deform}
\alias{predict.deform}
\title{Predict from a fitted \code{deform} object}
\usage{
\method{predict}{deform}(object, newdata = NULL, ...)
}
\arguments{
\item{object}{a fitted \code{deform} object}

\item{newdata}{a 2-column matrix of x and y coordinates}

\item{...}{currently just a placeholder}
}
\value{
A 2-column matrix of predicted x and y points for deformations and
a (2 + q)-column matrix for q-dimensional expansions.
}
\description{
Predict from a fitted \code{deform} object
}
\examples{

\donttest{

# fit a deformation model
data(solar)
m0 <- deform(solar$x, solar$z, solar$n)

# predict D-space points for original locations
predict(m0)

}

# predictions for one-dimensional expansion model with specified locations 
# and standard error estimates
data(solar)
m1 <- expand(solar$x, solar$z, solar$n)
xvals <- seq(-123.3, -122.2, by = .1)
yvals <- seq(49, 49.4, by = .1)
xyvals <- expand.grid(xvals, yvals)
predict(m1, xyvals, se.fit = TRUE)

}
