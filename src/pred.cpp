// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double maxtau = 0.9996;
const double twomaxtau = 2.0 * maxtau;

// [[Rcpp::export(.pcovfn)]]
arma::mat pcovfn(Rcpp::List pars, arma::mat z, int n){

double p1 = pars[0];
double p2 = pars[1];
double p3 = pars[2];
  
double sigma = exp(p1);
double tau = maxtau / (1.0 + exp(-p2));
double halfdelta = 1.0 / (1.0 + exp(-p3));
double sill = sigma * tau;
double dsq;

arma::mat C(n, n);
C.fill(R_NaN);

arma::vec dz;

for (int i=0; i < n; i++) {

  for (int j=i; j < n; j++) {

    dz = z.row(i).t() - z.row(j).t();
    dsq = sum(dz % dz);

    if (i == j) {
      C(i, j) = sigma;
    } else {
      if (dsq > 0) {
        C(i, j) = sill * exp(-R_pow(dsq, halfdelta));
        C(j, i) = C(i, j);
      }
    }
  }
}

return(C);

}

// [[Rcpp::export(.pcoscovfn)]]
arma::mat pcoscovfn(arma::vec pars, arma::mat z, int n){
  
double p0 = pars[0];
double p1 = pars[1];
double p2 = pars[2];
double p3 = pars[3];
  
arma::mat C(n, n);
C.fill(R_NaN);
  
double phi = exp(p0);
double sigma = exp(p1);
double tau = maxtau / (1.0 + exp(-p2));
double halfdelta = 1.0 / (1.0 + exp(-p3));
double sill = sigma * tau;
double dsq;

arma::vec dz;
          
for (int i=0; i < n; i++) {
            
  for (int j=i; j < n; j++) {
            
    dz = z.row(i).t() - z.row(j).t();
    dsq = sum(dz % dz);
      
    if (i == j) {
      C(i, j) = sigma;
    } else {
      if (dsq > 0) {
        C(i, j) = sill * exp(-R_pow(dsq, halfdelta)) * cos(sqrt(dsq) / phi);
        C(j, i) = C(i, j);
      }
    }
  }
}

return(C);
  
}

// [[Rcpp::export(.punitcoscovfn)]]
arma::mat punitcoscovfn(Rcpp::List pars, arma::mat z, int n){
  
double p0 = pars[0];
double p1 = 0.0;
double p2 = pars[1];
double p3 = pars[2];
  
arma::mat C(n, n);
C.fill(R_NaN);
  
double phi = exp(p0);
double sigma = exp(p1);
double tau = maxtau / (1.0 + exp(-p2));
double halfdelta = 1.0 / (1.0 + exp(-p3));
double sill = sigma * tau;
double dsq;
        
arma::vec dz;
        
for (int i=0; i < n; i++) {
          
  for (int j=i; j < n; j++) {
    dz = z.row(i).t() - z.row(j).t();
    dsq = sum(dz % dz);

    if (i == j) {
      C(i, j) = sigma;
    } else {
      if (dsq > 0) {
        C(i, j) = sill * exp(-R_pow(dsq, halfdelta)) * cos(sqrt(dsq) / phi);
        C(j, i) = C(i, j);
      }
    }
  }
}
        
return(C);
  
}

// [[Rcpp::export(.punitcovfn)]]
arma::mat punitcovfn(Rcpp::List pars, arma::mat z, int n){
  
double p2 = pars[0];
double p3 = pars[1];
  
arma::mat C(n, n);
C.fill(R_NaN);
  
double sigma = 1.0;
double tau = maxtau / (1.0 + exp(-p2));
double halfdelta = 1.0 / (1.0 + exp(-p3));
double sill = sigma * tau;
double dsq;
      
arma::vec dz;
      
for (int i=0; i < n; i++) {
        
  for (int j=i; j < n; j++) {
    dz = z.row(i).t() - z.row(j).t();
    dsq = sum(dz % dz);

    if (i == j) {
      C(i, j) = sigma;
    } else {
      if (dsq > 0) {
        C(i, j) = sill * exp(-R_pow(dsq, halfdelta));
        C(j, i) = C(i, j);
      }
    }
  }
}

return(C);
  
}

// [[Rcpp::export(.pdampedcoscovfn)]]
arma::mat pdampedcoscovfn(arma::vec pars, arma::mat z, int n){
  
double p0 = pars[0];
double p1 = pars[1];
double p2 = pars[2];
  
arma::mat C(n, n);
C.fill(R_NaN);
  
double phi = exp(p0);
double sigma = exp(p1);
double tau = maxtau / (1.0 + exp(-p2));
double sill = sigma * tau;
double dsq;

arma::vec dz;
          
for (int i=0; i < n; i++) {
            
  for (int j=i; j < n; j++) {
            
    dz = z.row(i).t() - z.row(j).t();
    dsq = sum(dz % dz);
      
    if (i == j) {
      C(i, j) = sigma;
    } else {
      if (dsq > 0) {
        C(i, j) = sill * exp(-sqrt(dsq)) * cos(sqrt(dsq) / phi);
        C(j, i) = C(i, j);
      }
    }
  }
}

return(C);
  
}

