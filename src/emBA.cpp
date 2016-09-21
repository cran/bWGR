#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emBA(NumericVector y, NumericMatrix gen,
          double df = 5, double R2 = 0.5, int it = 75){
  int p = gen.ncol();
  int n = gen.nrow();
  double ve = 1;
  double va = 1;
  double vy = var(y);
  NumericVector xx(p);
  NumericVector vx(p);
  for(int k=0; k<p; k++){
    xx[k] = sum(gen(_,k)*gen(_,k));
    vx[k] = var(gen(_,k));
  }
  NumericVector b(p);
  NumericVector vb = b+va;
  NumericVector Lmb = ve/vb;
  double MSx = sum(vx);
  double Sb = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double shape_hp  = 1.1;
  double rate_hp = (shape_hp-1)/Sb;
  double mu = mean(y);
  NumericVector e = y-mu;
  double b0,eM,h2,Sb_j;
  for(int i=0; i<it; i++){
    Sb_j = (1,p*df/2+shape_hp)/(sum(1/vb)/2+rate_hp);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e = e-gen(_,j)*(b[j]-b0);
      vb[j] = (Sb_j+b[j]*b[j]+va)/(df+2);
    }
    va = sum(b*b)/p;
    ve = (sum(e*e)+Se)/(n+df);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  double vg = sum(vb);
  h2 = vg/(vg+sqrt(ve));
  NumericVector fit(n);
  for(int k=0; k<n; k++){
    fit[k] = sum(gen(k,_)*b)+mu;
  }
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("h2") = h2);
}
