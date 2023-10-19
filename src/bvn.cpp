// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double logtwopi = log(2.0 * M_PI);
const double roottwopi = sqrt(2.0 * M_PI);

double dnrm_std(double x) {
  return(exp(-(x * x)/2) / roottwopi);
}

double dnrm(double x, double mu, double sig) {
  return(dnrm_std((x - mu)/sig)/sig);
}

double pnrm_std(double x) {
  return 0.5 * erfc(-x * M_SQRT1_2);
}

double pnrm(double x, double mu, double sig) {
  return pnrm_std((x - mu) / sig);
}

double pnorm2(double x) {
  double out = R::pnorm(x, 0.0, 1.0, 1, 0);
  return(out);}

double pbvn_neg(double DH, double DK, double R) {
  
  double TWOPI = 6.283185307179586;
  int I, J;
  double IS[2] = {-1, 1};
  
  int LG;
  double W[10], X[10];
  
  if (fabs(R) < .3) {
    LG = 3;
    W[0] = 0.1713244923791705;
    W[1] = 0.3607615730481384;
    W[2] = 0.4679139345726904;
    X[0] = -0.9324695142031522;
    X[1] = -0.6612093864662647;
    X[2] = -0.2386191860831970;
  }
  else if (fabs(R) < .75) {
    LG = 6;
    W[0] = .04717533638651177;
    W[1] = 0.1069393259953183;
    W[2] = 0.1600783285433464;
    W[3] = 0.2031674267230659;
    W[4] = 0.2334925365383547;
    W[5] = 0.2491470458134029;
    X[0] = -0.9815606342467191;
    X[1] = -0.9041172563704750;
    X[2] = -0.7699026741943050;
    X[3] = -0.5873179542866171;
    X[4] = -0.3678314989981802;
    X[5] = -0.1252334085114692;
  } 
  else {
    LG = 10;
    W[0] = .01761400713915212;
    W[1] = .04060142980038694;
    W[2] = .06267204833410906;
    W[3] = .08327674157670475;
    W[4] = 0.1019301198172404;
    W[5] = 0.1181945319615184;
    W[6] = 0.1316886384491766;
    W[7] = 0.1420961093183821;
    W[8] = 0.1491729864726037;
    W[9] = 0.1527533871307259;
    X[0] = -0.9931285991850949;
    X[1] = -0.9639719272779138;
    X[2] = -0.9122344282513259;
    X[3] = -0.8391169718222188;
    X[4] = -0.7463319064601508;
    X[5] = -0.6360536807265150;
    X[6] = -0.5108670019508271;
    X[7] = -0.3737060887154196;
    X[8] = -0.2277858511416451;
    X[9] = -0.07652652113349733;
  }
  
  
  double H = DH;
  double K = DK;
  double HK = H * K;
  double BVN = 0;
  
  double HS, ASR, SN, AS, A, B, BS, C, D, XS, RS;
  
  if (fabs(R) < .925) {
    if (fabs(R) > 0) {
      HS = ( H*H + K*K )/2;
      ASR = asin(R);
      for(I = 0; I < LG; I++) {
        for (J = 0; J < 2; J++) {
          SN = sin( ASR*(  IS[J]*X[I] + 1 )/2 );
          BVN = BVN + W[I]*exp( ( SN*HK-HS )/( 1-SN*SN ) );
        }
      }
      BVN = BVN*asin(R)/( 2*TWOPI );
    }
    //   BVN = BVN + pnorm(-H, 0.0, 1.0, 1, 0)*pnorm(-K, 0.0, 1.0, 1, 0);
    BVN = BVN + pnorm2(-H)*pnorm2(-K);
  } else {
    if (R < 0) {
      K = -K;
      HK = -HK;
    }
    if (fabs(R) < 1) {
      AS = ( 1 - R )*( 1 + R );
      A = sqrt(AS);
      BS = pow(H - K, 2);
      C = ( 4 - HK )/8;
      D = ( 12 - HK )/16;
      ASR = -( BS/AS + HK )/2;
      if (ASR > -100) BVN = A*exp(ASR)*( 1 - C*( BS - AS )*( 1 - D*BS/5 )/3 + C*D*AS*AS/5 );
      if (-HK < 100) {
        B = sqrt(BS);
        //       BVN = BVN - exp( -HK/2 )*sqrt(TWOPI)*pnorm(-B/A, 0.0, 1.0, 1, 0)*B*( 1 - C*BS*( 1 - D*BS/5 )/3 ); 
        BVN = BVN - exp( -HK/2 )*sqrt(TWOPI)*pnorm2(-B/A)*B*( 1 - C*BS*( 1 - D*BS/5 )/3 ); 
      }
      A = A/2;
      for(I = 0; I < LG; I++) {
        for (J = 0; J < 2; J++) {
          XS = pow(A*(  IS[J]*X[I] + 1 ), 2);
          RS = sqrt( 1 - XS );
          ASR = -( BS/XS + HK )/2;
          if (ASR > -100) {
            BVN = BVN + A*W[I]*exp( ASR )*( exp( -HK*( 1 - RS )/( 2*( 1 + RS ) ) )/RS - ( 1 + C*XS*( 1 + D*XS ) ) );
          }
        }
      }
      BVN = -BVN/TWOPI;
    }
    if (R > 0) {
      //     BVN =  BVN + pnorm(-fmax( H, K ), 0.0, 1.0, 1, 0);
      BVN =  BVN + pnorm2(-fmax( H, K ));
    } else {
      BVN = -BVN; 
      //     if (K > H) BVN = BVN + pnorm(K, 0.0, 1.0, 1, 0) - pnorm(H, 0.0, 1.0, 1, 0);
      if (K > H) BVN = BVN + pnorm2(K) - pnorm2(H);
    }
  }
  
  return(BVN);
  
}

double pbvn(double x, double y, double rho) {
  return pbvn_neg(-x, -y, rho);
}

double f00(double x, double y, double p0, double rho) {
  double sigma = exp(p0);
  return(-log(pbvn(x / sigma, y / sigma, rho)));
}

double f01(double x, double y, double p0, double rho) {
  double sigma = exp(p0);
  return(-log(pnrm(x, rho * y, sqrt(1 - rho * rho) * sigma)) - log(dnrm(x, 0, sigma)));
}

double f10(double x, double y, double p0, double rho) {
  return(f01(y, x, p0, rho));
}

double f11(double x, double y, double p0, double rho) {
  double temp = 1 - rho * rho;
  double z = x * x - 2 * rho * x * y + y * y;
  return(.5 * z / temp / exp(2 * p0) + logtwopi + 2 * p0 + .5 * log(temp));
}

double nllh_bvn_censored_ij(double x, double y, double ux, double uy, double p0, double rho) {
  
  double out = 0.0;
  
  if (x > ux) {
    if (y > uy) {
      out = f11(x, y, p0, rho);
    } else {
      out = f10(x, uy, p0, rho);
    }
  } else {
    if (y > uy) {
      out = f01(ux, y, p0, rho);
    } else {
      out = f00(ux, uy, p0, rho);
    }
  }
  return out;
}


// [[Rcpp::export(.nllh_bvn_censored_ogram)]]
double nllh_bvn_censored_ogram(arma::vec pars, arma::vec x, arma::vec y, arma::vec ux, arma::vec uy, int freq) {
  
  int n=x.size();
  double out = 0.0;
  
  for (int i=0; i < n; i++) {
    
    if (i % freq == 0) {
      
      if (!R_IsNA(x[i])) {
      
        if (!R_IsNA(y[i])) {
        
          out += nllh_bvn_censored_ij(x[i], y[i], ux[i], uy[i], pars[0], pars[1]);
         
        }
        
      }
      
    }
    
  }
  
  return(out);
  
}
