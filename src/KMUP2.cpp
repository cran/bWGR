#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP KMUP2(NumericMatrix X, NumericVector b,
            NumericVector xx, NumericVector E,
            NumericVector L, int p, double Ve,
            double pi, IntegerVector ro){
              
  RNGScope scope;
  NumericVector g = b;
  NumericVector e = E;
  NumericVector d(p);
  double G;
  double ee;
  double EE;
  double include;
  double exclude;
  int i;
  
  for(int j=0; j<p; j++){
    
      i = ro(j);
      G = g[i];
      g[i] = R::rnorm(
      (sum(X(_,i)*E) + xx(i)*G)/(xx(i)+L(i)),
      sqrt(Ve/(xx(i)+L(i))));
    
    if(pi>0){
      
      e = E - X(_,i) * (g[i]-G);
      ee = sum(e*e);
      EE = sum(E*E);
      include = (1-pi) * exp(-ee/(2*Ve));
      exclude = pi * exp(-EE/(2*Ve));
      d[i] = R::rbinom(1,include/(include+exclude));
      E = e;
      
    }else{
      
      d[i] = 1;
      E = E - X(_,i) * (g[i]-G);
      
    }
  }
  
   return List::create(Named("b") = g,
                       Named("d") = d,
                       Named("e") = E);
}
