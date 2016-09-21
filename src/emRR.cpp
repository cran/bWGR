#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emRR(NumericVector y, NumericMatrix gen,
           double df = 5, double R2 = 0.5, int it = 75){
  int p = gen.ncol();
  int n = gen.nrow();
  double Lmb = 1;
  double vb = 1;
  double va = 1;
  double ve = 1;
  double vy = var(y);
  NumericVector xx(p);
  NumericVector vx(p);
  for(int k=0; k<p; k++){
    xx[k] = sum(gen(_,k)*gen(_,k));
    vx[k] = var(gen(_,k));
  }
  double MSx = sum(vx);
  double Sb = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double mu = mean(y);
  NumericVector b(p);
  NumericVector e = y-mu;
  double b0,eM,h2;
  for(int i=0; i<it; i++){
    for(int j=0; j<p; j++){
      b0 = b[j];
      b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e = e-gen(_,j)*(b[j]-b0);
    }
    vb = (sum(b*b)+Sb)/(p+df+2);
    ve = (sum(e*e)+Se)/(n+df+2);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  va = sqrt(vb)*mean(xx);
  ve = sqrt(ve);
  h2 = (va/(va+ve))*(va/(va+ve));
  NumericVector fit(n);
  for(int k=0; k<n; k++){
    fit[k] = sum(gen(k,_)*b)+mu;
  }
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("h2") = h2);
}
