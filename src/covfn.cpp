// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double maxtau = 0.9996;
const double twomaxtau = 2.0 * maxtau;

// [[Rcpp::export(.covfn)]]
arma::mat covfn(Rcpp::List pars, Rcpp::List X, int n){

double p1 = Rcpp::as<double>(pars[0]);
double p2 = Rcpp::as<double>(pars[1]);
double p3 = Rcpp::as<double>(pars[2]);
    
arma::mat C(n, n);
C.fill(R_NaN);

  if (p1 < 100) {
    if (p2 > -100) {
      if (p3 > -100) {

double sigma = exp(p1);
double tau = maxtau / (1.0 + exp(-p2));
double halfdelta = 1.0 / (1.0 + exp(-p3));
double sill = sigma * tau;
double dsq = 0.0;

int nz = X.size();

arma::mat z(n, nz);
arma::vec dz;

if (nz > 0) {
    for (int l=0; l < nz; l++) {
    z.col(l) = Rcpp::as<arma::mat>(X[l]) * Rcpp::as<arma::vec>(pars[3 + l]);
}
}

for (int i=0; i < n; i++) {

for (int j=i; j < n; j++) {

if (nz > 0) {
    dz = z.row(i).t() - z.row(j).t();
    dsq = sum(dz % dz);
}
  
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

}}}

return(C);

}

// [[Rcpp::export(.d1covfn)]]
Rcpp::List d1covfn(Rcpp::List pars, Rcpp::List X, int n){
  
  double ee13 = exp(Rcpp::as<double>(pars[0]));
  double ee8 = exp(-Rcpp::as<double>(pars[1]));
  double ee1 = exp(-Rcpp::as<double>(pars[2]));
  
  double ee2 = 1.0 + ee1;
  double ee6 = 1.0 / ee2;
  double ee9 = 1.0 + ee8;
  double ee16 = ee9 * ee2;
  
  double ee18 = ee13 / ee16;
  double ee19 = ee13 / ee9;
  double ee20 = ee8 / ee9;
  double ee22 = ee13 * ee1 / (ee9 * ee2 * ee2);
  
  double dsq = 0.0;
  double ee7, ee10, ee12, ee17, ee21;
  
  int nz = X.size();
  
  arma::mat z(n, nz);
  arma::vec dz;
  
  if (nz > 0) {
    for (int l=0; l < nz; l++) {
      z.col(l) = Rcpp::as<arma::mat>(X[l]) * Rcpp::as<arma::vec>(pars[3 + l]);
    }
  }
  
  arma::mat dp1(n, n);
  arma::mat dp2(n, n);
  arma::mat dp3(n, n);
  arma::cube dzn(n, n, nz);
  
  for (int i=0; i < n; i++) {
    
    for (int j=i; j < n; j++) {
      
      if (i == j) {
        dp1(i, j) = ee13;
        dp2(i, j) = 0.0;
        dp3(i, j) = 0.0;
        for (int l=0; l < nz; l++) dzn(i, j, l) = 0.0;
        
      } else {
        
        if (nz > 0) {
          dz = z.row(j).t() - z.row(i).t();
          dsq = sum(dz % dz);
        }
        
        ee7 = dsq;
        ee10 = R_pow(ee7, ee6);
        ee12 = exp(-ee10);
        ee17 = R_pow(ee7, ee6 - 1.0);
        ee21 = ee17 * ee12;
        
        dp1(i, j) = maxtau * ee19 * ee12;
        dp2(i, j) = dp1(i, j) * ee20;
        dp3(i, j) = - maxtau * ee10 * ee12 * ee22 * log(dsq); 
        for (int l=0; l < nz; l++) dzn(i, j, l) = -(2.0 * maxtau * (dz[l] * ee21 * ee18));
        
        // the other triangle
        
        dp1(j, i) = dp1(i, j);
        dp2(j, i) = dp2(i, j);
        dp3(j, i) = dp3(i, j);
        for (int l=0; l < nz; l++) dzn(j, i, l) = -dzn(i, j, l);
        
      }
    }
  }
  
  Rcpp::List out(4 + nz);
  
  arma::mat dzmat;
  
  out[0] = dp1;
  out[1] = dp2;
  out[2] = dp3;
  for (int l=0; l < nz; l++) {
    dzmat = dzn.slice(l);
    out[3 + l] = dzmat;
  }
  
  return(out);
  
}

// [[Rcpp::export(.d2covfn)]]
Rcpp::List d2covfn(Rcpp::List pars, Rcpp::List X, int n){

double p1 = Rcpp::as<double>(pars[0]);
double p2 = Rcpp::as<double>(pars[1]);
double p3 = Rcpp::as<double>(pars[2]);

double ee1 = exp(-p3);
double ee2 = 1.0 + ee1;
double ee3 = 1.0;
double ee8 = 1.0/ee2;
double ee10 = exp(-p2);
double ee11 = 1.0 + ee10;
double ee12 = ee8 - 1.0;
double ee16 = exp(p1);
double ee17 = ee11 * ee2;
double ee19 = ee17;
double ee22 = ee2 * ee2;
double ee24 = ee19 * ee19;
double ee28 = ee11 * ee22;
double ee29 = ee11 * ee11;
double ee30 = 2.0/ee2;
double ee64 = ee28 * ee28;

double ff1 = ee11 * ee2 * ee2 * ee2;
double ff2 = ee11 / ee24;
double ff4 = - (ee1 / ee28);
double ff5 = ee10/ee29;
double ff7 = (2.0 * (ee10/ee11) - 1.0) * ee10/ee29;
double ff8 = -(ee22 * ee10 * ee1/ee64);
double ff9 = -(2.0 * (ee2 * ee10/ee24));
double ff10 = 2.0 * (ee17 * ee1/ee64);
double ff11 = ee1/ee22;

double ee7, ee9;
double ee13, ee15, ee18;
double ee21, ee23, ee27;
double ee32, ee34;
double ee47;
double ee56;
double ee78;
double ee81, ee83, ee84, ee85, ee86, ee87, ee88, ee89;
double gg1, gg2, gg3;
double dsq = 0.0;

int nz = X.size();

arma::mat z(n, nz);
arma::vec dz;

if (nz > 0) {
  for (int l=0; l < nz; l++) {
    z.col(l) = Rcpp::as<arma::mat>(X[l]) * Rcpp::as<arma::vec>(pars[3 + l]);
  }
}

arma::mat dp1p1(n, n);
arma::mat dp1p2(n, n);
arma::mat dp1p3(n, n);

arma::cube dp1z(n, n, nz);

arma::mat dp2p2(n, n);
arma::mat dp2p3(n, n);

arma::cube dp2z(n, n, nz);

arma::mat dp3p3(n, n);

arma::cube dp3z(n, n, nz);

arma::cube dzlzm(n, n, nz * nz);

arma::cube dzn(n, n, nz);

for (int i=0; i < n; i++) {
  
  for (int j=i; j < n; j++) {
    
    if (i == j) {
      
      dp1p1(j, i) = exp(p1);
      dp1p2(j, i) = 0.0;
      dp1p3(j, i) = 0.0;
      
      for (int l=0; l < nz; l++) dp1z(i, j, l) = 0.0;
      
      dp2p2(j, i) = 0.0;
      dp2p3(j, i) = 0.0;
      
      for (int l=0; l < nz; l++) dp2z(i, j, l) = 0.0;
      
      dp3p3(j, i) = 0.0;
      
      for (int l=0; l < nz; l++) dp3z(i, j, l) = 0.0;
      
      for (int l=0; l < nz; l++) {
        for (int m=l; m < nz; m++) {
          dzlzm(i, j, nz * l + m) = 0.0;
        }
      }
      
    } else {
      
      if (nz > 0) {
        dz = z.row(j).t() - z.row(i).t();
        dsq = sum(dz % dz);
      }
      
      ee7 = dsq;
      ee9 = ee7;
      ee13 = R_pow(ee9, ee8);
      ee15 = exp(-ee13);
      ee18 = R_pow(ee9, ee12);
      ee21 = log(ee7);
      ee23 = ee12 * R_pow(ee9, ee8 - 2.0);
      ee27 = R_pow(ee9, 2.0 * ee12)/ee2;
      ee32 = 2.0 * ee23 - 2.0 * ee27;
      ee34 = R_pow(ee9, ee30 - 1.0);
      ee47 = ff2 * ee18 + (ee18 - ee34) * ee21/ ff1;
      gg1 = ee15 * ee16;
      gg2 = ee18 * gg1;
      ee56 = ee13 * gg1;
      ee78 = ee56 * ff4 * ee21;
      ee81 = gg1 * ff5;
      
      dp1p1(i, j) = maxtau * gg1 / ee11;
      dp1p2(i, j) = maxtau * ee81; 
      dp1p3(i, j) = maxtau * ee78;
      ee83 = -2.0 * (gg2 / ee19);
      
      for (int l=0; l < nz; l++) dp1z(i, j, l) = dz[l] * ee83;
      
      dp2p2(i, j) = maxtau * ff7 * gg1;
      dp2p3(i, j) = maxtau * ee13 * gg1 * ee21 * ff8;
      ee84 = gg2 * ff9;
      
      for (int l=0; l < nz; l++) dp2z(i, j, l) = dz[l] * ee84;
      
      gg3 = ee1 * gg1;
      dp3p3(i, j) = -(maxtau * ((ff10 * ee13 + ((ee13 - R_pow(ee9, ee30)) * ee21 * ff11 - ee13)/ee28) * gg3 * ee21));
      ee85 = -(2.0 * (ee47 * gg3));
      
      for (int l=0; l < nz; l++) dp3z(i, j, l) = dz[l] * ee85;
      
      ee88 = gg1 / ee19;
      ee86 = -(2.0 * ee18 * ee88);
      ee89 = ee32 * ee88 / ee3;
      ee87 = -(2.0 * ee89);
      
      for (int l=0; l < nz; l++) {
        for (int m=l; m < nz; m++) {
          if (l == m) {
            dzlzm(i, j, nz * l + m) = maxtau * (ee86 + dz[l] * dz[l] * ee87);
          } else {
            dzlzm(i, j, nz * l + m) = -(2.0 * maxtau * (dz[l] * dz[m] * ee89));
          }
        }
      }
      
      // the other triangle
      
      dp1p1(j, i) = dp1p1(i, j);
      dp1p2(j, i) = dp1p2(i, j);
      dp1p3(j, i) = dp1p3(i, j);
      
      for (int l=0; l < nz; l++) dp1z(j, i, l) = -dp1z(i, j, l);
      
      dp2p2(j, i) = dp2p2(i, j);
      dp2p3(j, i) = dp2p3(i, j);
      
      for (int l=0; l < nz; l++) dp2z(j, i, l) = -dp2z(i, j, l);
      
      dp3p3(j, i) = dp3p3(i, j);
      
      for (int l=0; l < nz; l++) dp3z(j, i, l) = -dp3z(i, j, l);
      
      for (int l=0; l < nz; l++) {
        for (int m=l; m < nz; m++) {
          dzlzm(j, i, nz * l + m) = dzlzm(i, j, nz * l + m);
        }
      }
      
    }
  }
}

Rcpp::List dp1(3 + nz);
Rcpp::List dp2(3 + nz);
Rcpp::List dp3(3 + nz);
Rcpp::List dz1l(3 + nz);
Rcpp::List dz2l(3 + nz);

dp1[0] = dp1p1;
dp1[1] = dp1p2;
dp1[2] = dp1p3;

dp2[1] = dp2p2;
dp2[2] = dp2p3;

dp3[2] = dp3p3;

Rcpp::List dzl(nz);

for (int l=0; l < nz; l++) {
  Rcpp::List dzll(nz + 3);
  dp1[3 + l] = dp1z.slice(l);
  dp2[3 + l] = dp2z.slice(l);
  dp3[3 + l] = dp3z.slice(l);
  for (int m=l; m < nz; m++) dzll[3 + m] = dzlzm.slice(l * nz + m);
  dzl[l] = dzll;
}

Rcpp::List out(4 + nz);

out[0] = dp1;
out[1] = dp2;
out[2] = dp3;
for (int l=0; l < nz; l++) out[3 + l] = dzl[l];

return(out);

}
