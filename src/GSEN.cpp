#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP GSEN(NumericMatrix X, NumericVector b,
          NumericVector O,
          NumericVector xx, NumericVector e,
          double L, double A, int p){
  NumericVector g(p);
  NumericVector E = e;
  NumericVector OLS(p);
  NumericVector L1(p);
  NumericVector L2(p);
  double Lmb1 = L*(1-A);
  double Lmb2 = L*A;

  //for(int i=0; i<p; i++){
  int i;
  for(int j=0; j<p; j++){
    i = O(j)-1;
    
    if(A==1){
      g[i] = (sum(X(_,i)*E)+xx(i)*b[i])/(xx(i)+L);
    }else{
      OLS(i) = (sum(X(_,i)*E)+xx(i)*b[i])/(xx(i));
      if(OLS(i)>0){
        L1(i) = OLS(i) - (Lmb1/(2*xx(i)));
        if(L1(i)>0){
          L2(i) = xx(i)/(Lmb2+xx(i));
          g[i] = L1(i)*L2(i);
        }else{
          g[i] = 0;
        }
      }else{
        L1(i) = OLS(i) + (Lmb1/(2*xx(i)));
        if(L1(i)<0){
          L2(i) = xx(i)/(Lmb2+xx(i));
          g[i] = L1(i)*L2(i);
        }else{
          g[i] = 0;
        }
      }
      
    }
    E = E - X(_,i)*(g[i]-b[i]);
  }
  
  return List::create(Named("b")=g,Named("e")=E);
}
