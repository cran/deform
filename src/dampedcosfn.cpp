// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double maxtau = 0.9996;
const double twomaxtau = 2.0 * maxtau;

// [[Rcpp::export(.dampedcoscovfn)]]
arma::mat dampedcoscovfn(Rcpp::List pars, Rcpp::List X, int n){

double p0 = Rcpp::as<double>(pars[0]);
double p1 = Rcpp::as<double>(pars[1]);
double p2 = Rcpp::as<double>(pars[2]);

arma::mat C(n, n);
C.fill(R_NaN);

if (p0 < 100) {
  if (p1 < 100) {
    if (p2 > -100) {

double phi = exp(p0);
double sigma = exp(p1);
double tau = maxtau / (1.0 + exp(-p2));
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
      C(i, j) = sill * exp(-sqrt(dsq)) * cos(sqrt(dsq) / phi);
      C(j, i) = C(i, j);
    }
  }
}
}

}}}

return(C);

}

// [[Rcpp::export(.d1dampedcoscovfn)]]
Rcpp::List d1dampedcoscovfn(Rcpp::List pars, Rcpp::List X, int n){

double p0 = Rcpp::as<double>(pars[0]);
double p1 = Rcpp::as<double>(pars[1]);
double p2 = Rcpp::as<double>(pars[2]);
double dsq = 0.0;

// fixed stuff

double ee4 = exp(p0);
double ee7 = exp(-p2);
double ee9 = 1 + ee7;
double ee11 = exp(p1);

double ee3, ee5, ee10, ee12, ee13, ee14, ee15, ee16, eed;

int nz = X.size();

arma::mat z(n, nz);
arma::vec dz;

if (nz > 0) {
    for (int l=0; l < nz; l++) {
    z.col(l) = Rcpp::as<arma::mat>(X[l]) * Rcpp::as<arma::vec>(pars[3 + l]);
}
}

arma::mat dp0(n, n);
arma::mat dp1(n, n);
arma::mat dp2(n, n);
arma::cube dzn(n, n, nz);

for (int i=0; i < n; i++) {

for (int j=i; j < n; j++) {

if (i == j) {
  dp0(i, j) = 0.0;
  dp1(i, j) = ee11;
  dp2(i, j) = 0.0;
  for (int l=0; l < nz; l++) dzn(i, j, l) = 0.0;
  
  } else {
      
  if (nz > 0) {
      dz = z.row(j).t() - z.row(i).t();
      dsq = sum(dz % dz);
  }
  
  ee3 = sqrt(dsq);
  ee5 = ee3/ee4;
  ee10 = exp(-ee3);
  ee12 = cos(ee5);
  ee13 = sin(ee5);
  ee14 = ee9 * ee3;
  ee15 = ee12 + ee13/ee4;
  ee16 = maxtau * ee12;
  
dp0(i, j) = maxtau * ee10 * ee11 * ee13 * ee3/(ee9 * ee4);
dp1(i, j) = ee16 * ee10 * ee11/ee9;
dp2(i, j) = ee16 * ee7 * ee10 * ee11/ (ee9 * ee9);

eed = -(maxtau * ee15 * ee10 * ee11/ee14);

  for (int l=0; l < nz; l++) dzn(i, j, l) = dz[l] * eed;

// the other triangle

dp0(j, i) = dp0(i, j);
dp1(j, i) = dp1(i, j);
dp2(j, i) = dp2(i, j);
for (int l=0; l < nz; l++) dzn(j, i, l) = -dzn(i, j, l);
  
}
}
}

Rcpp::List out(3 + nz);

arma::mat dzmat;

out[0] = dp0;
out[1] = dp1;
out[2] = dp2;
for (int l=0; l < nz; l++) {
    dzmat = dzn.slice(l);
    out[3 + l] = dzmat;
}

return(out);

}

// [[Rcpp::export(.d2dampedcoscovfn)]]
Rcpp::List d2dampedcoscovfn(Rcpp::List pars, Rcpp::List X, int n){

double p0 = Rcpp::as<double>(pars[0]);
double p1 = Rcpp::as<double>(pars[1]);
double p2 = Rcpp::as<double>(pars[2]);
double dsq = 0.0;

double ee4 = exp(p0);
double ee7 = exp(-p2);
double ee10 = 1 + ee7;
double ee13 = exp(p1);
double ee15 = ee10 * ee10;
double ee16 = ee10 * ee4;

double ee2, ee3, ee5, ee8, ee9;
double ee12, ee14, ee17;
double ee20, ee21, ee22, ee23, ee27, ee28;
double ee39, ee43, ee49, ee54;
double ff1, ff2, ff3;
int nz = X.size();

arma::mat z(n, nz);
arma::vec dz;

if (nz > 0) {
    for (int l=0; l < nz; l++) {
    z.col(l) = Rcpp::as<arma::mat>(X[l]) * Rcpp::as<arma::vec>(pars[3 + l]);
}
}

arma::mat dp0p0(n, n);
arma::mat dp0p1(n, n);
arma::mat dp0p2(n, n);

arma::cube dp0z(n, n, nz);

arma::mat dp1p1(n, n);
arma::mat dp1p2(n, n);

arma::cube dp1z(n, n, nz);

arma::mat dp2p2(n, n);

arma::cube dp2z(n, n, nz);

arma::cube dzlzm(n, n, nz * nz);

arma::cube dzn(n, n, nz);

for (int i=0; i < n; i++) {

for (int j=i; j < n; j++) {

if (i == j) {

dp0p0(j, i) = 0.0;
dp0p1(j, i) = 0.0;
dp0p2(j, i) = 0.0;

for (int l=0; l < nz; l++) dp0z(i, j, l) = 0.0;

dp1p1(j, i) = ee13;
dp1p2(j, i) = 0.0;

for (int l=0; l < nz; l++) dp1z(i, j, l) = 0.0;

dp2p2(j, i) = 0.0;

for (int l=0; l < nz; l++) dp2z(i, j, l) = 0.0;

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

ee2 = dsq;
ee3 = sqrt(ee2);
ee5 = ee3/ee4;
ee8 = cos(ee5);
ee9 = sin(ee5);
ee12 = exp(-ee3);
ee14 = ee8 + ee9/ee4;
ee17 = ee8/ee4;
// ee18 = d1 * maxtau;
// ee19 = d2 * maxtau;
ee20 = ee16 * ee3;
ee21 = ee10 * ee3;
ee22 = ee15 * ee3;
ee23 = (ee17 - 2 * ee9)/ee4;
// ee24 = ee18 * ee14;
// ee25 = ee19 * ee14;
ee27 = ee8 * ee3/ee4;
ee28 = maxtau * ee8;
// ee29 = -(d1 * d2 * maxtau * (ee23 - (ee14/ee3 + ee8)) * ee12 * ee13/(ee10 * ee2));
// ee30 = -(ee24 * ee7 * ee12 * ee13/ee22);
// ee31 = -(ee24 * ee12 * ee13/ee21);
// ee32 = -(ee25 * ee7 * ee12 * ee13/ee22);
// ee33 = -(ee25 * ee12 * ee13/ee21);
// ee35 = (1 - ee3) * ee9 + ee27;
// ee37 = ee23 - ee8;
ee39 = (ee9 - ee17) * ee3 - ee9;
ee43 = ee28 * ee7 * ee12 * ee13/ee15;
ee49 = maxtau * ee7 * ee12 * ee13 * ee9 * ee3/(ee15 * ee4);
ee54 = maxtau * ee12 * ee13 * ee9 * ee3/ee16;
  
// p0

dp0p0(i, j) = -(maxtau * (ee27 + ee9) * ee12 * ee13 * ee3/ee16);
dp0p1(i, j) = ee54;
dp0p2(i, j) = ee49;

ff1 = -(maxtau * ee39 * ee12 * ee13/ee20);
for (int l=0; l < nz; l++) dp0z(i, j, l) = dz[l] * ff1;

// p1

dp1p1(i, j) = ee28 * ee12 * ee13/ee10;
dp1p2(i, j) = ee43;

ff2 = -(maxtau * ee14 * ee12 * ee13/ee21);
for (int l=0; l < nz; l++) dp1z(i, j, l) = dz[l] * ff2;

// p2

dp2p2(i, j) = maxtau * (2 * (ee7/ee10) - 1) * ee8 * ee7 * ee12 * ee13/ee15;

ff3 = -(maxtau * ee14 * ee7 * ee12 * ee13/ee22);
for (int l=0; l < nz; l++) dp2z(i, j, l) = dz[l] * ff3;

// d1, d2

for (int l=0; l < nz; l++) {
  for (int m=l; m < nz; m++) {
    if (l == m) {
      dzlzm(i, j, nz * l + m) = -(maxtau * (((ee17 - 2 * ee9)/ee4 - ee8) * R_pow(dz[l]/ee3, 2) + 0.5 * ((2 - 0.5 * (R_pow(2 * dz[l], 2)/ee2)) * ee14/ee3)) * ee12 * ee13/ee10);
    } else {
      dzlzm(i, j, nz * l + m) = -(dz[l] * dz[m] * maxtau * (ee23 - (ee14/ee3 + ee8)) * ee12 * ee13/(ee10 * ee2));
    }
  }
}

// the other triangle

dp0p0(j, i) = dp0p0(i, j);
dp0p1(j, i) = dp0p1(i, j);
dp0p2(j, i) = dp0p2(i, j);

for (int l=0; l < nz; l++) dp0z(j, i, l) = -dp0z(i, j, l);

dp1p1(j, i) = dp1p1(i, j);
dp1p2(j, i) = dp1p2(i, j);

for (int l=0; l < nz; l++) dp1z(j, i, l) = -dp1z(i, j, l);

dp2p2(j, i) = dp2p2(i, j);

for (int l=0; l < nz; l++) dp2z(j, i, l) = -dp2z(i, j, l);

for (int l=0; l < nz; l++) {
    for (int m=l; m < nz; m++) {
        dzlzm(j, i, nz * l + m) = dzlzm(i, j, nz * l + m);
    }
}

}
}
}

Rcpp::List dp0(3 + nz);
Rcpp::List dp1(3 + nz);
Rcpp::List dp2(3 + nz);

dp0[0] = dp0p0;
dp0[1] = dp0p1;
dp0[2] = dp0p2;

dp1[1] = dp1p1;
dp1[2] = dp1p2;

dp2[2] = dp2p2;

Rcpp::List dzl(nz);

for (int l=0; l < nz; l++) {
    Rcpp::List dzll(nz + 3);
    dp0[3 + l] = dp0z.slice(l);
    dp1[3 + l] = dp1z.slice(l);
    dp2[3 + l] = dp2z.slice(l);
    for (int m=l; m < nz; m++) dzll[3 + m] = dzlzm.slice(l * nz + m);
    dzl[l] = dzll;
}
  
Rcpp::List out(3 + nz);

out[0] = dp0;
out[1] = dp1;
out[2] = dp2;
for (int l=0; l < nz; l++) out[3 + l] = dzl[l];

return(out);

}
