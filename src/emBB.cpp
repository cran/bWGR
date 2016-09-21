#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP emBB(NumericVector y, NumericMatrix gen,
          double df = 5, double R2 = 0.5,
          int it = 75, double Pi = 0.95){
  int p = gen.ncol();
  int n = gen.nrow();
  double va = 1;
  double ve = 1;
  NumericVector d(p);
  NumericVector b(p);
  NumericVector vb = b+va;
  NumericVector Lmb = ve/vb;
  double vy = var(y);
  if(Pi>0.5){
    Pi = 1-Pi;
  } 
  NumericVector xx(p);
  NumericVector vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));
  }
  double MSx = sum(vx)*Pi;
  double Sb = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double shape_hp  = 1.1;
  double rate_hp = (shape_hp-1)/Sb;
  double mu = mean(y);
  NumericVector e = y-mu;
  NumericVector e1(n);
  NumericVector e2(n);
  double b0,b1,LR,eM,h2,Sb_j,C;
  double Pi0 = (1-Pi)/Pi;
  double MD = Pi;
  for(int i=0; i<it; i++){
    C = -0.5/ve;
    Sb_j = (1,p*df/2+shape_hp)/(sum(1/vb)/2+rate_hp);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e1 = e-gen(_,j)*(b1-b0);
      e2 = e-gen(_,j)*(0-b0);
      LR = Pi0*exp(C*(sum(e2*e2)-sum(e1*e1)));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j]/MD;
      vb[j] = (Sb_j+b[j]*b[j]+va)/(df+2);
      e = e - gen(_,j)*(b1-b0);
    }
    MD = max(d);
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

