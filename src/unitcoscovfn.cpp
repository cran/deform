// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double maxtau = 0.9996;
const double twomaxtau = 2.0 * maxtau;

// [[Rcpp::export(.unitcoscovfn)]]
arma::mat unitcoscovfn(Rcpp::List pars, Rcpp::List X, int n){

double p0 = Rcpp::as<double>(pars[0]);
double p1 = 0.0;
double p2 = Rcpp::as<double>(pars[1]);
double p3 = Rcpp::as<double>(pars[2]);
    
arma::mat C(n, n);
C.fill(R_NaN);

if (p0 < 100) {
  if (p2 > -100) {
      if (p3 > -100) {

double phi = exp(p0);
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
      C(i, j) = sill * exp(-R_pow(dsq, halfdelta)) * cos(sqrt(dsq) / phi);
      C(j, i) = C(i, j);
    }
  }
}
}

}}}

return(C);

}

// [[Rcpp::export(.d1unitcoscovfn)]]
Rcpp::List d1unitcoscovfn(Rcpp::List pars, Rcpp::List X, int n){

double p0 = Rcpp::as<double>(pars[0]);
double p1 = 0.0;
double p2 = Rcpp::as<double>(pars[1]);
double p3 = Rcpp::as<double>(pars[2]);
double dsq = 0.0;

double ee4 = exp(-p3);
double ee5 = 1 + ee4;
double ee6 = exp(p0);
double ee8 = 1/ee5;
double ee10 = ee5 * ee5;
double ee12 = exp(-p2); 
double ee14 = 1 + ee12; 
double ee16 = exp(p1);  

double ee2, ee7, ee9;
double ee11, ee15, ee17, ee18, ee19;
double ee22, ee23, ee24, ee25;

int nz = X.size();

arma::mat z(n, nz);
arma::vec dz;

if (nz > 0) {
    for (int l=0; l < nz; l++) {
    z.col(l) = Rcpp::as<arma::mat>(X[l]) * Rcpp::as<arma::vec>(pars[3 + l]);
}
}

arma::mat dp0(n, n);
// arma::mat dp1(n, n);
arma::mat dp2(n, n);
arma::mat dp3(n, n);
arma::cube dzn(n, n, nz);

for (int i=0; i < n; i++) {

for (int j=i; j < n; j++) {

if (i == j) {
  dp0(i, j) = 0.0;
  // dp1(i, j) = ee16;
  dp2(i, j) = 0.0;
  dp3(i, j) = 0.0;
  for (int l=0; l < nz; l++) dzn(i, j, l) = 0.0;
  
  } else {
      
  if (nz > 0) {
      dz = z.row(j).t() - z.row(i).t();
      dsq = sum(dz % dz);
  }

ee2 = dsq;
ee7 = sqrt(ee2);
ee9 = ee7/ee6;   
ee11 = R_pow(ee2, ee8);  
ee15 = exp(-ee11);   
ee17 = cos(ee9); 
ee18 = maxtau * ee17;
ee19 = sin(ee9); 
ee22 = 0.5 * (ee19/(ee6 * ee7)) + ee17 * R_pow(ee2, ee8 - 1)/ee5;  
ee23 = ee18 * ee15;
ee25 = ee16 / ee14;

dp0(i, j) = maxtau * ee15 * ee19 * ee7 * ee25 / ee6;
// dp1(i, j) = ee23 * ee25;
dp2(i, j) = ee23 * ee12 * ee25 / ee14;
dp3(i, j) = -(ee18 * ee11 * ee15 * ee4 * log(ee2) * ee25 / ee10);
if (Rcpp::traits::is_nan<REALSXP>(dp3(i, j))) {
  Rcpp::Rcout << "p3 is " << p3 << std::endl;
  Rcpp::Rcout << "ee18 is " << ee18 << std::endl;
  Rcpp::Rcout << "ee11 is " << ee11 << std::endl;
  Rcpp::Rcout << "ee15 is " << ee15 << std::endl;
  Rcpp::Rcout << "ee4 is " << ee4 << std::endl;
  Rcpp::Rcout << "ee25 is " << ee25 << std::endl;
  Rcpp::Rcout << "ee2 is " << ee2 << std::endl;
}
ee24 = -(2 * (maxtau * ee22 * ee15 * ee25));
for (int l=0; l < nz; l++) dzn(i, j, l) = dz[l] * ee24;

// the other triangle

dp0(j, i) = dp0(i, j);
// dp1(j, i) = dp1(i, j);
dp2(j, i) = dp2(i, j);
dp3(j, i) = dp3(i, j);
for (int l=0; l < nz; l++) dzn(j, i, l) = -dzn(i, j, l);
  
}
}
}

Rcpp::List out(3 + nz);

arma::mat dzmat;

out[0] = dp0;
// out[1] = dp1;
out[1] = dp2;
out[2] = dp3;
for (int l=0; l < nz; l++) {
    dzmat = dzn.slice(l);
    out[3 + l] = dzmat;
}

return(out);

}

// [[Rcpp::export(.d2unitcoscovfn)]]
Rcpp::List d2unitcoscovfn(Rcpp::List pars, Rcpp::List X, int n){

double p0 = Rcpp::as<double>(pars[0]);
double p1 = 0.0;
double p2 = Rcpp::as<double>(pars[1]);
double p3 = Rcpp::as<double>(pars[2]);
double dsq = 0.0;

double ee4 = exp(-p3);
double ee5 = exp(p0);
double ee6 = 1 + ee4;
double ee8 = 1/ee6;
double ee11 = exp(-p2);
double ee15 = ee8 - 1;
double ee16 = 1 + ee11;
double ee20 = exp(p1);
double ee26 = ee6 * ee6;
double ee27 = ee16 * ee16;
double ee28 = ee16 * ee26;
double ee34 = ee5 * ee5;
double ee35 = ee16 * ee5;
double ee45 = 2/ee6;

double ff1, ff3, ff4, ff5, ff6;
int nz = X.size();

arma::mat z(n, nz);
arma::vec dz;

if (nz > 0) {
    for (int l=0; l < nz; l++) {
    z.col(l) = Rcpp::as<arma::mat>(X[l]) * Rcpp::as<arma::vec>(pars[3 + l]);
}
}

arma::mat dp0p0(n, n);
// arma::mat dp0p1(n, n);
arma::mat dp0p2(n, n);
arma::mat dp0p3(n, n);

arma::cube dp0z(n, n, nz);

// arma::mat dp1p1(n, n);
// arma::mat dp1p2(n, n);
// arma::mat dp1p3(n, n);
// 
// arma::cube dp1z(n, n, nz);

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

dp0p0(j, i) = 0.0;
// dp0p1(j, i) = 0.0;
dp0p2(j, i) = 0.0;
dp0p3(j, i) = 0.0;

for (int l=0; l < nz; l++) dp0z(i, j, l) = 0.0;

// dp1p1(j, i) = ee20;
// dp1p2(j, i) = 0.0;
// dp1p3(j, i) = 0.0;
// 
// for (int l=0; l < nz; l++) dp1z(i, j, l) = 0.0;

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
  
double ee2, ee7, ee9, ee12, ee13, ee14, ee18, ee19;
double ee21, ee22, ee24, ee25, ee29;
double ee32, ee33, ee37;
double ee43, ee44, ee52, ee53, ee58;
double ee65, ee66, ee74;
  
ee2 = dsq;
ee7 = sqrt(ee2);
ee9 = ee7/ee5;
ee12 = R_pow(ee2, ee8);
ee13 = cos(ee9);
ee14 = sin(ee9);
ee18 = R_pow(ee2, ee15);
ee19 = exp(-ee12);
ee21 = ee5 * ee7;
ee22 = ee13 * ee18;
ee24 = 0.5 * (ee14/ee21) + ee22/ee6;
ee25 = log(ee2);
ee29 = ee18 * ee14;
ee32 = maxtau * ee13;
ee33 = ee21 * ee21;
ee37 = (ee15 * ee13 * R_pow(ee2, ee8 - 2) - (ee24 * ee18 + 0.5 * 
  (ee29/ee21)))/ee6 + 0.5 * (0.5 * (ee13/(ee2 * ee34)) - 
  0.5 * (ee5 * ee14/(ee33 * ee7)));
ee43 = ee32 * ee12 * ee19;
ee44 = maxtau * ee19;
// ee46 = ee32 * ee19;
ee52 = -(maxtau * ee12 * ee19 * ee4 * ee20 * ee25 * ee14 * 
  ee7/(ee28 * ee5));
ee53 = -(ee43 * ee11 * ee4 * ee20 * ee25/(ee27 * ee26));
// ee54 = -(ee43 * ee4 * ee20 * ee25/ee28);
// ee56 = ((ee18 - R_pow(ee2, ee45 - 1)) * ee13/ee6 - 0.5 * (ee12 * 
// ee14/ee21)) * ee25 + ee22;
ee58 = (ee18 + ee18 * ee25/ee6) * ee13 - ee24 * ee12 * ee25;
// ee60 = (0.5 * (ee13/ee21) - ee29/ee6) * ee7 + 0.5 * (ee14/ee7);
ee65 = ee29 * ee7/(ee6 * ee34) - 0.5 * (ee13/R_pow(ee5, 3) + ee14 * 
  ee7/ee33);
ee66 = 2 * ee24;
// ee69 = ee46 * ee11 * ee20/ee27;
ee74 = ee44 * ee11 * ee20 * ee14 * ee7/(ee27 * ee5);
// ee78 = ee44 * ee20 * ee14 * ee7/ee35;

dp0p0(i, j) = -(maxtau * (ee13 * ee7/ee5 + ee14) * ee19 * 
        ee20 * ee7/ee35);
// dp0p1(i, j) = ee78;
dp0p2(i, j) = ee74;
dp0p3(i, j) = ee52;

ff1 = -(2 * (maxtau * ee65 * ee19 * ee5 * ee20/ee16));
for (int l=0; l < nz; l++) dp0z(i, j, l) = dz[l] * ff1;

// dp1p1(i, j) = ee46 * ee20/ee16;
// dp1p2(i, j) = ee69;
// dp1p3(i, j) = ee54;
// 
// ff2 = -(2 * (maxtau * ee24 * ee19 * ee20/ee16));
// for (int l=0; l < nz; l++) dp1z(i, j, l) = dz[l] * ff2;

dp2p2(i, j) = maxtau * (2 * (ee11/ee16) - 1) * ee13 * ee19 * ee11 * ee20/ee27;
dp2p3(i, j) = ee53;

ff3 = -(2 * (maxtau * ee24 * ee19 * ee11 * ee20/ee27));
for (int l=0; l < nz; l++) dp2z(i, j, l) = dz[l] * ff3;

dp3p3(i, j) = -(maxtau * ((ee12 - R_pow(ee2, ee45)) * R_pow(ee4/ee26, 2) * 
                ee25 + (2 * (ee4/ee6) - 1) * ee12 * ee4/ee26) * 
                ee13 * ee19 * ee20 * ee25/ee16);

ff4 = -(2 * (maxtau * ee58 * ee19 * ee4 * ee20/ee28));
for (int l=0; l < nz; l++) dp3z(i, j, l) = dz[l] * ff4;

ff5 = -(maxtau * ee19 * ee20/ee16);
ff6 = -(4 * (maxtau * ee37 * ee19 * ee20/ee16));
for (int l=0; l < nz; l++) {
    for (int m=l; m < nz; m++) {
        if (l == m) {
            dzlzm(i, j, nz * l + m) = ff5 * (ee37 * 4 * dz[l] * dz[l] + ee66);
        } else {
            dzlzm(i, j, nz * l + m) = dz[l] * dz[m] * ff6;
        }
    }
}

// the other triangle

dp0p0(j, i) = dp0p0(i, j);
// dp0p1(j, i) = dp0p1(i, j);
dp0p2(j, i) = dp0p2(i, j);
dp0p3(j, i) = dp0p3(i, j);

for (int l=0; l < nz; l++) dp0z(j, i, l) = -dp0z(i, j, l);

// dp1p1(j, i) = dp1p1(i, j);
// dp1p2(j, i) = dp1p2(i, j);
// dp1p3(j, i) = dp1p3(i, j);
// 
// for (int l=0; l < nz; l++) dp1z(j, i, l) = -dp1z(i, j, l);

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

Rcpp::List dp0(3 + nz);
// Rcpp::List dp1(3 + nz);
Rcpp::List dp2(3 + nz);
Rcpp::List dp3(3 + nz);

dp0[0] = dp0p0;
// dp0[1] = dp0p1;
dp0[1] = dp0p2;
dp0[2] = dp0p3;

// dp1[1] = dp1p1;
// dp1[2] = dp1p2;
// dp1[3] = dp1p3;

dp2[1] = dp2p2;
dp2[2] = dp2p3;

dp3[2] = dp3p3;

Rcpp::List dzl(nz);

for (int l=0; l < nz; l++) {
    Rcpp::List dzll(nz + 3);
    dp0[3 + l] = dp0z.slice(l);
    // dp1[4 + l] = dp1z.slice(l);
    dp2[3 + l] = dp2z.slice(l);
    dp3[3 + l] = dp3z.slice(l);
    for (int m=l; m < nz; m++) dzll[3 + m] = dzlzm.slice(l * nz + m);
    dzl[l] = dzll;
}
  
Rcpp::List out(3 + nz);

out[0] = dp0;
// out[1] = dp1;
out[1] = dp2;
out[2] = dp3;
for (int l=0; l < nz; l++) out[3 + l] = dzl[l];

return(out);

}
