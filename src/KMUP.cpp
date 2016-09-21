#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP KMUP(NumericMatrix X, NumericVector b,
          NumericVector d, NumericVector xx,
          NumericVector E, NumericVector L,
          double Ve, double pi){
  
  RNGScope scope;
  int p = X.ncol();
  NumericVector g = b;
  NumericVector e1 = E;
  NumericVector e2 = E;
  double G,G0,Gp,D,pj,LR;
  double Cons = -0.5/Ve;
  double Pi0 = pi/(1-pi);
  
  for(int i=0; i<p; i++){
    G0 = g[i];
    G = R::rnorm((sum(X(_,i)*E)+G0)/(xx(i)+L(i)),sqrt(Ve/(xx(i)+L(i))));
    Gp = G*pi;
    e1 = E - X(_,i)*(G-G0);
    if(pi>0){
      e2 = E - X(_,i)*(Gp-G0);
      LR = Pi0*exp(Cons*(sum(e1*e1)-sum(e2*e2)));
      pj = 1-1/(1+LR);
      D = R::rbinom(1,sqrt(pj));
     if(D==0){
      d[i] = 0;
      g[i] = Gp;
      E = e2;
     }else{
      d[i] = 1;
      g[i] = G;
      E = e1;
     }
    }else{
      d[i] = 1;
      g[i] = G;
      E = e1;
    }
  }
   return List::create(Named("b") = g,
                       Named("d") = d,
                       Named("e") = E);
}
