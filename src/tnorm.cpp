// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double gaussconst = 0.5 * log(2.0 * arma::datum::pi);

// [[Rcpp::export(.tnormr)]]
double tnormr(arma::vec pars, arma::vec yvec, arma::vec rvec)
{
    
double mu = pars[0];
double lsigma = pars[1];
double sigma = exp(lsigma);
int nobs = yvec.size();
double y, z, r, v;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
r = rvec[j];
z = (y - mu) / sigma;
v = (r - mu) / sigma;

nllh += 0.5 * z * z;
nllh += log(arma::normcdf(v));

}

nllh += nobs * (lsigma + gaussconst);

if (!arma::is_finite(nllh)) {
  nllh = 1.0e12;
}

return(nllh);

}

// [[Rcpp::export(.tnorml)]]
double tnorml(arma::vec pars, arma::vec yvec, arma::vec lvec)
{
  
  double mu = pars[0];
  double lsigma = pars[1];
  double sigma = exp(lsigma);
  int nobs = yvec.size();
  double y, z, l, w;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    l = lvec[j];
    z = (y - mu) / sigma;
    w = (l - mu) / sigma;
    
    nllh += 0.5 * z * z;
    nllh += log(1.0 - arma::normcdf(w));
    
  }
  
  nllh += nobs * (lsigma + gaussconst);
  
  if (!arma::is_finite(nllh)) {
    nllh = 1.0e12;
  }
  
  return(nllh);
    
}

// [[Rcpp::export(.tnormlr)]]
double tnormlr(arma::vec pars, arma::vec yvec, arma::vec lvec, arma::vec rvec)
{
  
  double mu = pars[0];
  double lsigma = pars[1];
  double sigma = exp(lsigma);
  int nobs = yvec.size();
  double y, z, r, l, v, w;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    l = lvec[j];
    r = rvec[j];
    z = (y - mu) / sigma;
    v = (l - mu) / sigma;
    w = (r - mu) / sigma;
    
    nllh += 0.5 * z * z;
    nllh += log(arma::normcdf(w) - arma::normcdf(v));
    
  }
  
  nllh += nobs * (lsigma + gaussconst);
  
  if (!arma::is_finite(nllh)) {
    nllh = 1.0e12;
  }
  
  return(nllh);
    
}
