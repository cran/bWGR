#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP KMUP2(NumericMatrix X, NumericVector Use, NumericVector b,
          NumericVector d, NumericVector xx, NumericVector E,
          NumericVector L, double Ve, double pi){
  
  RNGScope scope;
  int p = X.ncol();
  int n0 = X.nrow();
  int n = Use.size();
  double bg = n0/n;
  NumericVector E0(n);
  NumericVector H(n);
  for(int k=0; k<n; k++){
    E0[k] = E[Use[k]];
  }
  NumericVector g = b;
  NumericVector e1 = E0;
  NumericVector e2 = E0;
  double G,G0,Gp,D,pj,LR;
  double Cons = -0.5/Ve;
  double Pi0 = (1-pi)/pi;
  
  for(int i=0; i<p; i++){
    G0 = g[i];
    for(int x=0; x<n; x++){
      H[x] = X(Use[x],i);
    }
    G = R::rnorm((sum(H*E0)+G0)/(xx(i)*bg+L(i)),
                  sqrt(Ve/(xx(i)*bg+L(i))));
    Gp = G*pi;
    e1 = E0 - H*(G-G0);
    if(pi>0){
      e2 = E0 - H*(Gp-G0);
      LR = Pi0*exp(Cons*(sum(e1*e1)-sum(e2*e2)));
      pj = 1-1/(1+LR);
      D = R::rbinom(1,pj);
     if(D==0){
      d[i] = 0;
      g[i] = Gp;
      E0 = e2;
     }else{
      d[i] = 1;
      g[i] = G;
      E0 = e1;
     }
    }else{
      d[i] = 1;
      g[i] = G;
      E0 = e1;
    }
  }
  
  
   return List::create(Named("b") = g,
                       Named("d") = d,
                       Named("e") = E0);
}
