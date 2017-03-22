#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP KMUP(NumericMatrix X, NumericVector b, NumericVector d, NumericVector xx, NumericVector E, NumericVector L, double Ve, double pi){  
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
    G = R::rnorm((sum(X(_,i)*E)+G0)/(xx(i)+L(i)),
                 sqrt(Ve/(xx(i)+L(i))));
    Gp = G*pi;
    e1 = E - X(_,i)*(G-G0);
    if(pi>0){
      e2 = E - X(_,i)*(Gp-G0);
      LR = Pi0*exp(Cons*(sum(e1*e1)-sum(e2*e2)));
      pj = 1-1/(1+LR);
      D = R::rbinom(1,pj);
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

// [[Rcpp::export]]
SEXP KMUP2(NumericMatrix X, NumericVector Use, NumericVector b,  NumericVector d, NumericVector xx, NumericVector E, NumericVector L, double Ve, double pi){  
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
  double Pi0 = pi/(1-pi);
  
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

// [[Rcpp::export]]
SEXP emBA(NumericVector y, NumericMatrix gen, double df = 4, double R2 = 0.5){
  int it = 200;
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
    Sb_j = (p*df/2+shape_hp)/(sum(1/vb)/2+rate_hp);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e = e-gen(_,j)*(b[j]-b0);
      vb[j] = (Sb_j+b[j]*b[j]+va)/(df+2)+ve/(xx[j]+Lmb[j]);
    }
    va = sum(b*b)/p;
    ve = (sum(e*e)+Se)/(n+df);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  double vg = sum(vb);
  h2 = vg/(vg+ve);
  NumericVector fit(n);
  for(int k=0; k<n; k++){
    fit[k] = sum(gen(k,_)*b)+mu;
  }
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Va") = vg,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBB(NumericVector y, NumericMatrix gen, double df = 4, double R2 = 0.5, double Pi = 0.7){
  R2 = R2/5;
  int it = 200;
  df = df*10;
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
    Sb_j = (p*df/2+shape_hp)/(sum(1/vb)/2+rate_hp);
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
  double vg = sum(vb)/Pi;
  h2 = vg/(vg+ve);
  NumericVector fit(n);
  for(int k=0; k<n; k++){
    fit[k] = sum(gen(k,_)*b)+mu;
  }
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Va") = vg,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBC(NumericVector y, NumericMatrix gen, double df = 4, double R2 = 0.5, double Pi = 0.7){
  R2 = R2/5;
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  double Lmb = 1;
  double vb = 1;
  double va = 1;
  double ve = 1;
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
  double mu = mean(y);
  NumericVector b(p);
  NumericVector d(p);
  NumericVector e = y-mu;
  NumericVector e1(n);
  NumericVector e2(n);
  double b0,b1,LR,eM,h2,C;
  double Pi0 = (1-Pi)/Pi;
  double MD = Pi;
  for(int i=0; i<it; i++){
    C = -0.5/ve;
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e1 = e-gen(_,j)*(b1-b0);
      e2 = e-gen(_,j)*(0-b0);
      LR = Pi0*exp(C*(sum(e2*e2)-sum(e1*e1)));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j];
      e = e - gen(_,j)*(b1-b0);
    }
    MD = max(d);
    b = b/MD;
    vb = (sum(b*b)+Sb)/(p+df);
    ve = (sum(e*e)+Se)/(n+df);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  va = mean(xx)*(sum(b*b)+Sb)/(p+df)/Pi;
  h2 = va/(va+ve);
  NumericVector fit(n);
  for(int k=0; k<n; k++){
    fit[k] = sum(gen(k,_)*b)+mu;
  }
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Va") = va,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emRR(NumericVector y, NumericMatrix gen, double df = 4, double R2 = 0.5){
  int it = 200;
  df = df*10000;
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
                      Named("Va") = va,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBL(NumericVector y, NumericMatrix gen, double R2 = 0.5, double alpha = 0.02){
  int it = 200;
  double h2 = R2;  
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  NumericVector b(p);
  double b0,eM,Half_L2;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  // Regulation coefficients
  double cxx = mean(xx);
  double Lmb1 = cxx*((1-h2)/h2)*alpha*0.5;
  double Lmb2 = cxx*((1-h2)/h2)*(1-alpha);
  double OLS, G;
  // Loop
  for(int i=0; i<it; i++){
    // Regression coefficients loop
    for(int j=0; j<p; j++){
      // Ordinary Least Square
      b0 = b[j];
      OLS = (sum(gen(_,j)*e)+xx[j]*b0);
      Half_L2 = 0.5*OLS/(xx[j]+cxx);
      // Regularization of positive OLS
      if(OLS>0){
        G = 0.5*(OLS-Lmb1)/(Lmb2+xx(j));
        if(G>0){b[j] = G+Half_L2;}else{b[j] = Half_L2;}
      }else{
        // Regularization of negative OLS
        G = 0.5*(OLS+Lmb1)/(Lmb2+xx(j));
        if(G<0){b[j] = G+Half_L2;}else{b[j] = Half_L2;}
      }
      // Residuals update
      e = e-gen(_,j)*(b[j]-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;}
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  h2 = 1 - var(e)/var(y);
  // Output
  return List::create(Named("mu") = mu, Named("b") = b, Named("hat") = fit, Named("h2") = h2); }

// [[Rcpp::export]]
SEXP emDE(NumericVector y, NumericMatrix gen, double R2 = 0.5){
  // Convergence criteria
  int maxit = 300;
  double tol = 10e-6;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  double b0,eM;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  NumericVector vx(p);
  for(int k=0; k<p; k++){vx[k] = var(gen(_,k));}
  double cxx = sum(vx)*(1-R2)/R2;
  // Regulation coefficients
  double Ve;
  NumericVector Vb(p);
  NumericVector b(p);
  NumericVector Lmb = p+cxx;
  double b1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      // Ordinary Least Square
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(Lmb(j)+xx(j));
      b[j] = b1;
      // Residuals update
      e = e-gen(_,j)*(b1-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components
    Ve = sum(e*y)/(n-1);
    Vb = b*b+(Ve/(xx+Lmb));
    Lmb = sqrt(cxx*Ve/Vb);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("Vb")=Vb,
                      Named("Ve")=Ve,
                      Named("h2")=sum(Vb)/(sum(Vb)+Ve));
}

// [[Rcpp::export]]
SEXP emEN(NumericVector y, NumericMatrix gen, double R2 = 0.5, double alpha = 0.02){
  // Convergence criteria
  int maxit = 300;
  double tol = 10e-11;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  NumericVector b(p);
  double b0,eM;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  NumericVector vx(p);
  for(int k=0; k<p; k++){vx[k] = var(gen(_,k));}
  double cxx = sum(vx)*(1-R2)/R2;
  // Regulation coefficients
  double Ve, Va;
  double Sy = sd(y);
  double Lmb = cxx;
  double Lmb1 = 0.5*Lmb*alpha*Sy;
  double Lmb2 = Lmb*(1-alpha);
  double trAC22 = sum(1/(xx+Lmb));
  double OLS, b1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      // Ordinary Least Square
      b0 = b[j];
      OLS = (sum(gen(_,j)*e)+xx[j]*b0);
      // Regularization of positive OLS
      if(OLS>0){
        // Regularization of positive OLS
        b1 = (OLS-Lmb1)/(Lmb2+xx(j));if(b1<0){b1=0;};
      }else{
        // Regularization of negative OLS
        b1 = (OLS+Lmb1)/(Lmb2+xx(j));if(b1>0){b1=0;};
      }
      b[j] = b1;
      // Residuals update
      e = e-gen(_,j)*(b1-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components
    Ve = sum(e*y)/(n-1);
    Va = (sum(b*b)+trAC22*Ve)/p;
    Lmb = Ve/Va;
    Lmb1 = 0.5*Lmb*alpha*Sy;
    Lmb2 = Lmb*(1-alpha);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("Va")=Va*cxx,
                      Named("Ve")=Ve,
                      Named("h2")=Va*cxx/(Va*cxx+Ve));
}

