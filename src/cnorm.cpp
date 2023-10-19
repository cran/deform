// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double gaussconst = 0.5 * log(2.0 * arma::datum::pi);

// [[Rcpp::export(.cnorml)]]
double cnorml(arma::vec pars, arma::vec yvec, arma::vec lvec)
{
    
  double mu = pars[0];
  double lsigma = pars[1];
  double sigma = exp(lsigma);
  int nobs = yvec.size();
  double y, z, l, v;
  double nllh = 0.0;

  for (int j=0; j < nobs; j++) {

    y = yvec[j];
  
    if (arma::is_finite(y)) {
    
      l = lvec[j];
      z = (y - mu) / sigma;
      v = (l - mu) / sigma;

      if (y > l) {
  
        nllh += 0.5 * z * z + lsigma + gaussconst;

      } else {
  
        nllh -= log(arma::normcdf(v));

      }
  
    }
  
  }

  if (!arma::is_finite(nllh)) {
    nllh = 1.0e12;
  }

  return(nllh);

}

// [[Rcpp::export(.cnormr)]]
double cnormr(arma::vec pars, arma::vec yvec, arma::vec rvec)
{
  
  double mu = pars[0];
  double lsigma = pars[1];
  double sigma = exp(lsigma);
  int nobs = yvec.size();
  double y, z, r, w;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    
    if (arma::is_finite(y)) {
    
      r = rvec[j];
      z = (y - mu) / sigma;
      w = (r - mu) / sigma;
    
      if (y < r) {
      
        nllh += 0.5 * z * z + lsigma + gaussconst;
      
      } else {
      
        nllh -= log(1 - arma::normcdf(w));
    
      }
    }
  }

  if (!arma::is_finite(nllh)) {
    nllh = 1.0e12;
  }
  
  return(nllh);
    
}

// [[Rcpp::export(.cnormlr)]]
double cnormlr(arma::vec pars, arma::vec yvec, arma::vec lvec, arma::vec rvec)
{
  
  double mu = pars[0];
  double lsigma = pars[1];
  double sigma = exp(lsigma);
  int nobs = yvec.size();
  double y, z, r, l, v, w;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    
    if (arma::is_finite(y)) {
      
      l = lvec[j];
      r = rvec[j];
      z = (y - mu) / sigma;
      v = (l - mu) / sigma;
      w = (r - mu) / sigma;
    
      if (y < r) {
      
        if (y > l) {
      
          nllh += 0.5 * z * z + lsigma + gaussconst;
      
        } else {
       
          nllh -= log(arma::normcdf(v));
      
        }
      
      } else {
      
        nllh -= log(1 - arma::normcdf(w));
      
      }
      
    }  

  }
  
  if (!arma::is_finite(nllh)) {
    nllh = 1.0e12;
  }
  
  return(nllh);
    
}
