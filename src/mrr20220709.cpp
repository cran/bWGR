// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>

// [[Rcpp::export]]
SEXP mrr(Eigen::MatrixXd Y, Eigen::MatrixXd X){
  
  // Basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  int maxit = 1000;
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXd y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Sum of squares of X
  Eigen::MatrixXd XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXd XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx = XSX.colwise().sum();
  Eigen::VectorXd TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  Eigen::VectorXd ve = vy * 0.5;
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSx.array()).matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  
  // Beta tilde;
  Eigen::MatrixXd tilde = X.transpose() * y;
  
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Convergence control
  double deflate = 1.0, deflateMax = 0.75;
  Eigen::MatrixXd beta0(p,k), A(k,k);
  double cnv = 10.0, logtol = -10.0;
  int numit = 0;
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      // Update coefficient
      b0 = b.row(J)*1.0;
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      b.row(J) = b1;
      // Update residuals
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN.array();
    iVe = ve.array().inverse();
    
    // Genetic variance
    TildeHat = b.transpose()*tilde;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){ vb(i,i) = TildeHat(i,i)/TrXSX(i); }else{
        vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrXSX(i)+TrXSX(j));}}}
    
    // Bending
    A = vb*1.0;
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i!=j){A(i,j) *= deflate;}}}
    while(A.llt().info()==Eigen::NumericalIssue && deflate>deflateMax){
      Rcpp::Rcout << "Bend at "<< numit << "\n";  deflate = deflate - 0.01;
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*deflate;}}}
    }
    iG = A.inverse();
    
    // Print status
    cnv = log10((beta0.array()-b.array()).square().sum());  ++numit;
    if( numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; }
    if( numit % 200 == 0){ deflate = deflate - 0.01; }
    if( cnv<logtol ){break;}
    if(std::isnan(cnv)){ break;}
    
  }
  
  // Fitting the model
  h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXd hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Genetic correlations
  Eigen::MatrixXd GC(k,k);
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("Vb")=vb,
                            Rcpp::Named("Ve")=ve,
                            Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("bend")=deflate,
                            Rcpp::Named("cnv")=cnv);
  
}

// [[Rcpp::export]]
SEXP mrr_float(Eigen::MatrixXf Y, Eigen::MatrixXf X){
  
  // Basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  int maxit = 500;
  
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXf mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXf y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Sum of squares of X
  Eigen::MatrixXf XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXf XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXf MSx = XSX.colwise().sum();
  Eigen::VectorXf TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  Eigen::VectorXf ve = vy * 0.5;
  Eigen::VectorXf iVe = ve.array().inverse();
  Eigen::MatrixXf vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSx.array()).matrix().asDiagonal();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  
  // Beta tilde;
  Eigen::MatrixXf tilde = X.transpose() * y;
  
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Inner RGS
  std::vector<int> InnerRGSvec(k);
  for(int j=0; j<k; j++){InnerRGSvec[j]=j;}
  std::random_device rd2;
  std::mt19937 g2(rd2());
  int ri;
  
  // Convergence control
  Eigen::MatrixXf beta0(p,k);
  float cnv = 10.0, logtol = -8.0;
  int numit = 0;
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    std::shuffle(InnerRGSvec.begin(), InnerRGSvec.end(), g2);
    
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      b0 = b.row(J)*1.0;
      // Stranden and Garrick 2009
      LHS = vb * (XX.row(J).transpose().array() * iVe.array()).matrix().asDiagonal(); 
      for(int i=0; i<k; i++){ LHS(i,i) += 1.0; }
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = (vb * (RHS.array() * iVe.array()).matrix()).array();
      // Inner GS
      b1 = b.row(J)*1.0;
      for(int i=0; i<k; i++){
        ri = InnerRGSvec[i];
        b1(ri) = (RHS(ri)-(LHS.col(ri).array()*b1.array()).sum()+LHS(ri,ri)*b1(ri))/LHS(ri,ri);}
      // Update residuals
      b.row(J) = b1;
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN.array();
    iVe = ve.array().inverse();
    
    // Genetic variance
    TildeHat = b.transpose()*tilde;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){ vb(i,i) = TildeHat(i,i)/TrXSX(i); }else{
        vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrXSX(i)+TrXSX(j));}}}
    
    // Print status
    cnv = log10((beta0.array()-b.array()).square().colwise().sum().maxCoeff());
    ++numit; if( cnv<logtol ){break;}
    if(std::isnan(cnv)){ break;}
    
  }
  
  // Fitting the model
  h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXf hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Genetic correlations
  Eigen::MatrixXf GC(k,k);
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("Vb")=vb,
                            Rcpp::Named("Ve")=ve,
                            Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("cnv")=cnv);
  
}


// [[Rcpp::export]]
SEXP mrr2X(Eigen::MatrixXd Y, Eigen::MatrixXd X1, Eigen::MatrixXd X2){
  
  // Basic info
  int maxit = 1000;
  int k = Y.cols(), n0 = Y.rows(), p1 = X1.cols(), p2 = X2.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXd y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Sum of squares of X1
  Eigen::MatrixXd XX1(p1,k);
  for(int i=0; i<p1; i++){
    XX1.row(i) = X1.col(i).array().square().matrix().transpose() * Z;}
  
  // Sum of squares of X2
  Eigen::MatrixXd XX2(p2,k);
  for(int i=0; i<p2; i++){
    XX2.row(i) = X2.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX)1;
  Eigen::MatrixXd XSX1(p1,k);
  for(int i=0; i<p1; i++){
    XSX1.row(i) = XX1.row(i).transpose().array()*iN.array() - 
      ((X1.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx1 = XSX1.colwise().sum();
  Eigen::VectorXd TrXSX1 = n.array()*MSx1.array();
  
  // Compute Tr(XSX)2;
  Eigen::MatrixXd XSX2(p2,k);
  for(int i=0; i<p2; i++){
    XSX2.row(i) = XX2.row(i).transpose().array()*iN.array() - 
      ((X2.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx2 = XSX2.colwise().sum();
  Eigen::VectorXd TrXSX2 = n.array()*MSx2.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  // VE
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy = vy.array()*iN.array();
  Eigen::VectorXd ve = vy * 0.5;
  Eigen::VectorXd iVe = ve.array().inverse();
  // VB1
  Eigen::MatrixXd vb1(k,k), TildeHat1(k,k);
  vb1 = (ve.array()/MSx1.array()).matrix().asDiagonal();
  Eigen::MatrixXd iG1 = vb1.inverse();
  // VB2
  Eigen::MatrixXd vb2(k,k), TildeHat2(k,k);
  vb2 = (ve.array()/MSx2.array()).matrix().asDiagonal();
  Eigen::MatrixXd iG2 = vb2.inverse();
  
  // Beta tilde;
  Eigen::MatrixXd tilde1 = X1.transpose() * y;
  Eigen::MatrixXd tilde2 = X2.transpose() * y;
  
  // Initialize coefficient matrices
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd bA = Eigen::MatrixXd::Zero(p1,k);
  Eigen::MatrixXd bB = Eigen::MatrixXd::Zero(p2,k);
  
  // RGS
  std::vector<int> RGSvec1(p1);
  std::vector<int> RGSvec2(p2);
  for(int j=0; j<p1; j++){RGSvec1[j]=j;}
  for(int j=0; j<p2; j++){RGSvec2[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Convergence control
  double deflate1 = 1.0, deflate2 = 1.0, deflateMax = 0.75;
  Eigen::MatrixXd beta01(p1,k), beta02(p2,k), A(k,k);
  double cnv = 10.0, logtol = -10.0;
  int numit = 0;
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta01 = bA*1.0;
    beta02 = bB*1.0;
    
    // Randomized Gauss-Seidel loop 1
    std::shuffle(RGSvec1.begin(), RGSvec1.end(), g);
    for(int j=0; j<p1; j++){
      J = RGSvec1[j];
      // Update coefficient
      b0 = bA.row(J)*1.0;
      LHS = iG1;
      LHS.diagonal() += (XX1.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (X1.col(J).transpose()*e).array() + XX1.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      bA.row(J) = b1;
      // Update residuals
      e = (e-(X1.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Randomized Gauss-Seidel loop 2
    std::shuffle(RGSvec2.begin(), RGSvec2.end(), g);
    for(int j=0; j<p2; j++){
      J = RGSvec2[j];
      // Update coefficient
      b0 = bB.row(J)*1.0;
      LHS = iG2;
      LHS.diagonal() += (XX2.row(2).transpose().array() * iVe.array()).matrix();
      RHS = (X2.col(J).transpose()*e).array() + XX2.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      bB.row(J) = b1;
      // Update residuals
      e = (e-(X2.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN.array();
    iVe = ve.array().inverse();
    
    // Genetic variance 1
    TildeHat1 = bA.transpose()*tilde1;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){ vb1(i,i) = TildeHat1(i,i)/TrXSX1(i); }else{
        vb1(i,j) = (TildeHat1(i,j)+TildeHat1(j,i))/(TrXSX1(i)+TrXSX1(j));}}}
    // Bending 1
    A = vb1*1.0;
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i!=j){A(i,j) *= deflate1;}}}
    while(A.llt().info()==Eigen::NumericalIssue && deflate1>deflateMax){
      Rcpp::Rcout << "Bend V1 at "<< numit << "\n";
      deflate1 = deflate1 - 0.01;
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){ A(i,j) = vb1(i,j)*deflate1; }}}}
    iG1 = A.inverse();
    
    // Genetic variance 1
    TildeHat2 = bB.transpose()*tilde2;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){ vb2(i,i) = TildeHat2(i,i)/TrXSX2(i); }else{
        vb2(i,j) = (TildeHat2(i,j)+TildeHat2(j,i))/(TrXSX2(i)+TrXSX2(j));}}}
    // Bending 1
    A = vb2*1.0;
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i!=j){A(i,j) *= deflate2;}}}
    while(A.llt().info()==Eigen::NumericalIssue && deflate2>deflateMax){
      Rcpp::Rcout << "Bend V2 at "<< numit << "\n";
      deflate2 = deflate2 - 0.01;
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){ A(i,j) = vb2(i,j)*deflate2; }}}}
    iG2 = A.inverse();
    
    
    // Print status
    cnv = log10((beta01.array()-bA.array()).square().sum()) + log10((beta02.array()-bB.array()).square().sum());
    ++numit;
    if( numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){break;}
    if(std::isnan(cnv)){ break;}
    
  }
  
  // Fitting the model
  Eigen::MatrixXd hat = X1 * bA + X2 * bB;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Genetic correlations and genetic variance 1
  Eigen::MatrixXd GC1(k,k), va1(k,k);
  va1.diagonal() = vb1.diagonal().array() * MSx1.array();
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ GC1(i,j)=vb1(i,j)/(sqrt(vb1(i,i)*vb1(j,j)));}}
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ if(i!=j){ va1(i,j)= GC1(i,j)*sqrt(va1(i,i)*va1(j,j));}}}
  
  // Genetic correlations and genetic variance 2
  Eigen::MatrixXd GC2(k,k), va2(k,k);
  va2.diagonal() = vb2.diagonal().array() * MSx2.array();
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ GC2(i,j)=vb2(i,j)/(sqrt(vb2(i,i)*vb2(j,j)));}}
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ if(i!=j){ va2(i,j)= GC2(i,j)*sqrt(va2(i,i)*va2(j,j));}}}
  
  // Heritability
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  Eigen::VectorXd h2A = va1.diagonal().array() / ( va1.diagonal().array() + ve.array()).array();
  Eigen::VectorXd h2B = va2.diagonal().array() / ( va2.diagonal().array() + ve.array()).array();
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b1")=bA, Rcpp::Named("b2")=bB,
                            Rcpp::Named("hat")=hat, Rcpp::Named("h2")=h2,
                            Rcpp::Named("h2_1")=h2A, Rcpp::Named("h2_2")=h2B,
                            Rcpp::Named("GC1")=GC1, Rcpp::Named("GC2")=GC2,
                            Rcpp::Named("VE")=ve,
                            Rcpp::Named("VB1")=vb1,Rcpp::Named("VB2")=vb2,
                            Rcpp::Named("VA1")=va1, Rcpp::Named("VA2")=va2,
                            Rcpp::Named("bend1")=deflate1,Rcpp::Named("bend2")=deflate2,
                            Rcpp::Named("cnv")=cnv);}

// [[Rcpp::export]]
SEXP mrr_svd(Eigen::MatrixXd Y, Eigen::MatrixXd W){
  
  // Start setup
  int maxit = 500;
  double tol = 10e-9;
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), m = W.cols();
  
  // Center X
  Rcpp::Rcout << "Centering marker score matrix\n";
  Eigen::VectorXd xx = W.colwise().mean();
  for(int i=0; i<m; i++){ W.col(i) = W.col(i).array() - xx(i);}
  
  // Single value decomposition
  Rcpp::Rcout << "SVD of marker scores\n";
  Eigen::BDCSVD<Eigen::MatrixXd> svd(W, Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXd V = svd.matrixV();
  Eigen::MatrixXd X = svd.matrixU() * svd.singularValues().array().matrix().asDiagonal();
  int p = X.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  Rcpp::Rcout << "Centering Y\n";
  Eigen::VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXd y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();
  }
  
  // Sum of squares of X
  Rcpp::Rcout << "Computing diagonal elements of Z'Z\n";
  Eigen::MatrixXd XX(p,k);
  for(int i=0; i<p; i++){ XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  // Compute Tr(XSX);
  Eigen::MatrixXd XSX(p,k);
  for(int i=0; i<p; i++){ XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx = XSX.colwise().sum();
  Eigen::VectorXd TrXSX = n.array()*MSx.array();
  
  
  Rcpp::Rcout << "Set starting values for coefficients and variances\n";
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  Eigen::VectorXd ve = vy * 0.5;
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSx.array()).matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  // Beta tilde;
  Eigen::MatrixXd tilde = X.transpose() * y;
  Eigen::VectorXd TrDinvXSX(k);
  Eigen::MatrixXd Dinv(p,k);
  for(int i=0; i<k; i++){ XSX.col(i) = XSX.col(i).array() * n(i); }
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  // Bending and convergence control
  Eigen::MatrixXd A = vb*1.0, GC(k,k);
  double deflate = 1.0, deflateMax = 0.75;
  Eigen::MatrixXd beta0(p,k), vb0(k,k);
  Eigen::VectorXd CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  double cnv = 10.0;
  int numit = 0;
  double logtol = log10(tol);
  
  // Loop
  Rcpp::Rcout << "Starting Gauss-Seidel\n";
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    vb0 = vb*1.0;
    ve0 = ve*1.0;
    h20 = h2*1.0;
    
    // Randomized Gauss-Seidel loop
    for(int J=0; J<p; J++){
      // Update coefficient
      b0 = b.row(J)*1.0;
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() * iVe.array();
      b1 = LHS.llt().solve(RHS);
      b.row(J) = b1;
      // Update residuals
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN.array();
    iVe = ve.array().inverse();
    h2 = 1 - ve.array()/vy.array();
    
    // Get tilde-hat
    for(int i=0; i<k; i++){
        Dinv.col(i) = (XSX.col(i).array()/ve(i) + iG(i,i)).inverse().array();
        TrDinvXSX(i)  = (XSX.col(i).transpose() * Dinv.col(i));}
    TildeHat = b.transpose()* Dinv.cwiseProduct(tilde);
    
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ // Variances
          vb(i,i) = TildeHat(i,i)/TrDinvXSX(i);
        }else{ // Covariances
          vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrDinvXSX(i)+TrDinvXSX(j));
        }}}
    
    // Bending
    A = vb*1.0;
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i!=j){A(i,j) *= deflate;}}}
    while(A.llt().info()==Eigen::NumericalIssue && deflate>deflateMax){
      Rcpp::Rcout << "Bend at iteration "<< numit << "\n";
      deflate = deflate - 0.01;
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*deflate;}}}
    }    
    iG = A.inverse();
    
    // Covariances
    cnv = log10((beta0.array()-b.array()).square().sum());  CNV1(numit) = cnv;
    CNV2(numit) = log10((h20.array()-h2.array()).square().sum());
    CNV3(numit) = log10((vb0.array()-vb.array()).square().sum());
    
    // Print
    ++numit;
    if( numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){ Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n"; break; }
    if(std::isnan(cnv)){ break;}
    
  }
  
  Rcpp::Rcout << "Fitting final model\n";
  // Fitting the model
  Eigen::MatrixXd hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  Eigen::MatrixXd beta = V*b;
  
  // Correlations
  Rcpp::Rcout << "Estimating correlations\n";
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Resize convergence vectors
  Rcpp::Rcout << "Convergence statistics\n";
  Eigen::VectorXd CNV1b(numit),CNV2b(numit),CNV3b(numit);
  for(int i=0; i<numit; i++){ CNV1b(i)=CNV1(i);CNV2b(i)=CNV2(i);CNV3b(i)=CNV3(i);}
  
  // Null model Output
  Rcpp::List NullModelOutput = Rcpp::List::create(Rcpp::Named("Intercepts")=mu,
                                                  Rcpp::Named("MarkerEffects")=beta,
                                                  Rcpp::Named("FittedValues")=hat,
                                                  Rcpp::Named("Heritability")=h2,
                                                  Rcpp::Named("WCorrelations")=GC,
                                                  Rcpp::Named("VarBeta")=vb,
                                                  Rcpp::Named("VarResiduals")=ve,
                                                  Rcpp::Named("CovarianceBending")=deflate,
                                                  Rcpp::Named("ConvergenceBeta")=CNV1b,
                                                  Rcpp::Named("ConvergenceH2")=CNV2b,
                                                  Rcpp::Named("ConvergenceVar")=CNV3b,
                                                  Rcpp::Named("NumOfIterations")=numit);
  NullModelOutput.attr("class") = "WModel";
  return NullModelOutput;
  
}

// [[Rcpp::export]]
SEXP MRR3(Eigen::MatrixXd Y,
          Eigen::MatrixXd X,
          int maxit = 500,
          double tol = 10e-9,
          int cores = 1,
          bool TH = false,
          double NLfactor = 0.0,
          bool InnerGS = false,
          bool NoInv = false,
          bool HCS = false,
          bool XFA2 = false,
          int NumXFA = 2,
          double R2 = 0.5,
          double gc0 = 0.5, 
          double df0 = 0.0, 
          double weight_prior_h2 = 0.0,
          double weight_prior_gc = 0.0,
          double PenCor = 0.0,
          double MinCor = 1.0,
          double uncorH2below = 0.0,
          double roundGCupFrom = 1.0,
          double roundGCupTo = 1.0,
          double roundGCdownFrom = 1.0,
          double roundGCdownTo = 0.0,
          double bucketGCfrom = 1.0,
          double bucketGCto = 1.0,
          double DeflateMax = 0.9,
          double DeflateBy = 0.005,
          bool OneVarB = false,
          bool OneVarE = false,
          bool verbose = false){
  
  //Set multi-core processing
  if(cores!=1) Eigen::setNbThreads(cores);
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXd y(n0,k);
  for(int i=0; i<k; i++){y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Center X
  Eigen::VectorXd xx = X.colwise().mean();
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - xx(i);}
  
  // Sum of squares of X
  Eigen::MatrixXd XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXd XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx = XSX.colwise().sum();
  Eigen::VectorXd TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  
  Eigen::VectorXd ve = vy * (1-R2);
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  Eigen::VectorXd vbInit = ((vy*R2).array()/MSx.array());
  Eigen::VectorXd veInit = ve*1.0;
  vb = vbInit.array().matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  
  // Starting covariance values
  double tmp;
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      if(i>j){
        tmp = gc0 * sqrt(vb(i,i)*vb(j,j));
        vb(i,j) = tmp;
        vb(j,i) = tmp;
      }
    }
  }
  
  // Beta tilde;
  Eigen::MatrixXd tilde = X.transpose() * y;
  Eigen::VectorXd TrDinvXSX(k);
  Eigen::MatrixXd Dinv(p,k);
  if(TH){
    for(int i=0; i<k; i++){
      XSX.col(i) = XSX.col(i).array() * n(i);
    }
  }
  
  // Prior shape
  Eigen::MatrixXd Sb = vb*df0;
  Eigen::VectorXd Se = ve*df0;
  Eigen::VectorXd iNp = (n.array()+df0-1).inverse();
  
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  
  // Bending and convergence control
  Eigen::MatrixXd A = vb*1.0, GC(k,k);
  double bucketMean = 0.5*(bucketGCfrom+bucketGCto);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EVDofA(A);
  double MinDVb, inflate = 0.0, Deflate = 1.0;
  Eigen::MatrixXd beta0(p,k), vb0(k,k);
  Eigen::VectorXd CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  double cnv = 10.0;
  int numit = 0;
  double logtol = log10(tol);
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Inner RGS
  std::vector<int> InnerRGSvec(k);
  for(int j=0; j<k; j++){InnerRGSvec[j]=j;}
  std::random_device rd2;
  std::mt19937 g2(rd2());
  int ri;
  
  // Non-Linear weights for marker effects
  bool NonLinear = NLfactor!=0.0;
  Eigen::MatrixXd W(p,k);
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){  W(i,j) = 1.0; }}
  Eigen::VectorXd iVeWj = iVe*1.0;
  Eigen::VectorXd tmpW(p);
  double maxW, minW;
  
  // Objects for other variance structures
  double gs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(vb);
  Eigen::MatrixXd UDU(k,k);
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    vb0 = vb*1.0;
    ve0 = ve*1.0;
    h20 = h2*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    std::shuffle(InnerRGSvec.begin(), InnerRGSvec.end(), g2);
    
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      
      // System of equations - Traditional vs Stranden and Garrick 2009
      if(NoInv){
        LHS = vb * (XX.row(J).transpose().array() * iVeWj.array()).matrix().asDiagonal(); 
        for(int i=0; i<k; i++){ LHS(i,i) += 1.0; }
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = (vb * (RHS.array() * iVeWj.array()).matrix()).array();
      }else{
        LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVeWj.array()).matrix();
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = RHS.array() * iVeWj.array();
      }
      
      // Update coefficient
      b0 = b.row(J)*1.0;
      for(int i=0; i<k; i++){ iVeWj(i) = iVe(i)*W(J,i); }
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVeWj.array()   ).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() * iVeWj.array();
      
      // Inner GS
      if(InnerGS){
        b1 = b.row(J)*1.0;
        for(int i=0; i<k; i++){
          ri = InnerRGSvec[i];
          b1(ri) = (RHS(ri)-(LHS.col(ri).array()*b1.array()).sum()+LHS(ri,ri)*b1(ri))/LHS(ri,ri);}
      }else{
        b1 = LHS.llt().solve(RHS); 
      }
      
      // Update residuals
      b.row(J) = b1;
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Update marker weights
    if(NonLinear){
      W = b.cwiseAbs();
      for(int j=0; j<k; j++){
        maxW = W.col(j).maxCoeff(); minW = W.col(j).minCoeff();
        tmpW = NLfactor * (W.col(j).array()-minW)/(maxW-minW) + (1.0-NLfactor);
        tmpW = tmpW.array() + (1.0-tmpW.mean());
        W.col(j) = tmpW.array();
      }
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = (ve.array()+Se.array()) * iNp.array();
    h2 = 1 - ve.array()/vy.array();
    // Proportion-based prior
    if(weight_prior_h2>0){for(int i=0; i<k; i++){gs = ve(i)*(1-weight_prior_h2) + weight_prior_h2*veInit(i); ve(i) = gs*1.0;}}
    // Single variance
    if(OneVarE){tmp = ve.array().mean(); for(int i=0; i<k; i++) ve(i) = tmp*1.0;}
    iVe = ve.array().inverse();
    iVeWj = iVe*1.0;
    
    //Genetic variance
    
    // Get tilde-hat
    if(TH){
      for(int i=0; i<k; i++){
        Dinv.col(i) = (XSX.col(i).array()/ve(i) + iG(i,i)).inverse().array();
        TrDinvXSX(i)  = (XSX.col(i).transpose() * Dinv.col(i));
      }
      TildeHat = b.transpose()* Dinv.cwiseProduct(tilde);
    }else{
      TildeHat = b.transpose()*tilde;
    }
    
    // Estimate variances and covariance components
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ // Variances
          if(TH){
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrDinvXSX(i)+df0);
          }else{
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrXSX(i)+df0);
          }
        }else{ // Covariances
          if(TH){
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrDinvXSX(i)+TrDinvXSX(j)+df0);
          }else{
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrXSX(i)+TrXSX(j)+df0);
          }
        }}}
    
    
    if(weight_prior_h2>0){ // Proportion-based prior H2
      for(int i=0; i<k; i++){gs = vb(i,i)*(1-weight_prior_h2) + weight_prior_h2*vbInit(i); vb(i,i) = gs*1.0;}}
    if(weight_prior_gc>0){ // Proportion-based prior GC
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){
        if(i!=j){ GC(i,j) = (1.0-weight_prior_gc)*vb(i,j)/(sqrt(vb(i,i)*vb(j,j))) + gc0*weight_prior_gc;}else{GC(i,j) = 1.0;}}}
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){ if(i!=j){ vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}else{
        // Once calculation of GC without prior
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}}
      
      // Heterogeneous Compound Symmetry
      if(HCS){
        gs = 0.0;
        for(int i=0; i<k; i++){
          for(int j=0; j<k; j++){
            if(i>j){gs += GC(i,j);}}}
        gs = gs/((k*(k-1))/2);
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){ 
          if(i!=j){ GC(i,j) =  gs*1.0;}else{ GC(i,j) = 1.0; }}}
        // Extended Factor Analytics
      }else if(XFA2){
        es.compute(GC);
        UDU = es.eigenvalues()[k] * es.eigenvectors().col(k) * es.eigenvectors().col(k).transpose();
        for(int i=1; i<NumXFA; i++) UDU += es.eigenvalues()[k-i] * es.eigenvectors().col(k-i) * es.eigenvectors().col(k-i).transpose();
        GC = UDU * 1.0; for(int i=0; i<k; i++){ GC(i,i)=1.0; };
      }
      
      // Monkeying with the correlations
      for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
          if(i!=j){
            // Zero'ing  Correlations
            if(MinCor<1.0){ if(GC(i,j)<MinCor){ GC(i,j) = 0.0; }}
            // Penalize Correlations
            if(PenCor>0.0){  GC(i,j) = tanh(PenCor*abs(GC(i,j)))*GC(i,j);} 
            // Round Down
            if(roundGCdownFrom<1.0){ if(GC(i,j)<roundGCdownFrom){ GC(i,j) = roundGCdownTo*1.0; }}
            // Round Up
            if(roundGCupFrom<1.0){ if(GC(i,j)>roundGCupFrom){ GC(i,j) = roundGCupTo*1.0; }}
            // Bucket round
            if(bucketGCfrom<1.0){ if(GC(i,j)>bucketGCfrom && GC(i,j)<bucketGCto  ){ GC(i,j) =  bucketMean*1.0; }}
            // Min H2
            if(uncorH2below>0.0){ if(h2(i)<uncorH2below || h2(j)<uncorH2below  ){ GC(i,j) = 0.0; }}
          }}}
      
      // RECONSTRUCT COVARIANCE HERE AND ONLY ONCE
      // Single variance of beta
      if(OneVarB){tmp = TildeHat.diagonal().mean(); vb = GC * tmp;  }else{
        // Regular covariance reconstruction
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){
          if(i!=j){ vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}
      
      // Bending
      if( !NoInv || TH ){
        A = vb*1.0;
        for(int i=0; i<k; i++){ for(int j=0; j<k; j++){if(i!=j){A(i,j) *= Deflate;}}}
        if(A.llt().info()==Eigen::NumericalIssue && Deflate>DeflateMax){
          Deflate -= DeflateBy; if(verbose) Rcpp::Rcout <<  Deflate;
          for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*Deflate;}}}}
        EVDofA.compute(A);
        MinDVb = EVDofA.eigenvalues().minCoeff();
        if( MinDVb < 0.0 ){ 
          if(verbose) Rcpp::Rcout << ".";
          inflate = inflate - MinDVb*1.001;
          A.diagonal().array() += inflate;}
        iG = A.inverse();
      }
      
      // Compute convergence and print status
      
      //cnv = log10((beta0.array()-b.array()).square().sum());
      cnv = log10((beta0.array()-b.array()).square().colwise().sum().maxCoeff());
      CNV1(numit) = cnv; if(std::isnan(cnv)){ if(verbose){Rcpp::Rcout << "Numerical issue! Job aborted (it=" << numit << ")\n";} break;}
      CNV2(numit) = log10((h20.array()-h2.array()).square().sum());
      CNV3(numit) = log10((vb0.array()-vb.array()).square().sum());
      
      // Print
      ++numit;
      if( verbose && numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
      if( cnv<logtol ){ if(verbose){Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n";} break;}
      if( numit == maxit && verbose){ Rcpp::Rcout << "Model did not converge\n"; }    
    
  }
  
  // Fitting the model
  Eigen::MatrixXd hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Resize convergence vectors
  Eigen::VectorXd CNV1b(numit),CNV2b(numit),CNV3b(numit);
  for(int i=0; i<numit; i++){
    CNV1b(i)=CNV1(i);
    CNV2b(i)=CNV2(i);
    CNV3b(i)=CNV3(i);
  }
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("vb")=vb,
                            Rcpp::Named("ve")=ve,
                            Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("bend")=Deflate,
                            Rcpp::Named("cnvB")=CNV1b,
                            Rcpp::Named("cnvH2")=CNV2b,
                            Rcpp::Named("cnvV")=CNV3b,
                            Rcpp::Named("b_Weights")=W,
                            Rcpp::Named("Its")=numit);
}

// [[Rcpp::export]]
SEXP MRR3F(Eigen::MatrixXf Y,
          Eigen::MatrixXf X,
          int maxit = 500,
          float tol = 10e-9,
          int cores = 1,
          bool TH = false,
          float NonLinearFactor = 0.0,
          bool InnerGS = false,
          bool NoInversion = false,
          bool HCS = false,
          bool XFA = false,
          int NumXFA = 2,
          float prior_R2 = 0.5,
          float gc_prior_df = 0.5, 
          float var_prior_df = 0.0, 
          float weight_prior_h2 = 0.0,
          float weight_prior_gc = 0.0,
          float PenCor = 0.0,
          float MinCor = 1.0,
          float uncorH2below = 0.0,
          float roundGCupFrom = 1.0,
          float roundGCupTo = 1.0,
          float roundGCdownFrom = 1.0,
          float roundGCdownTo = 0.0,
          float bucketGCfrom = 1.0,
          float bucketGCto = 1.0,
          float DeflateMax = 0.9,
          float DeflateBy = 0.005,
          bool OneVarB = false,
          bool OneVarE = false,
          bool verbose = false){
  
  //Set multi-core processing
  if(cores!=1) Eigen::setNbThreads(cores);
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXf mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXf y(n0,k);
  for(int i=0; i<k; i++){y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Center X
  Eigen::VectorXf xx = X.colwise().mean();
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - xx(i);}
  
  // Sum of squares of X
  Eigen::MatrixXf XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXf XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXf MSx = XSX.colwise().sum();
  Eigen::VectorXf TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  
  Eigen::VectorXf ve = vy * (1-prior_R2);
  Eigen::VectorXf iVe = ve.array().inverse();
  Eigen::MatrixXf vb(k,k), TildeHat(k,k);
  Eigen::VectorXf vbInit = ((vy*prior_R2).array()/MSx.array());
  Eigen::VectorXf veInit = ve*1.0;
  vb = vbInit.array().matrix().asDiagonal();
  Eigen::MatrixXf iG = vb.inverse();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  
  // Starting covariance values
  float tmp;
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      if(i>j){
        tmp = gc_prior_df * sqrt(vb(i,i)*vb(j,j));
        vb(i,j) = tmp;
        vb(j,i) = tmp;
      }
    }
  }
  
  // Beta tilde;
  Eigen::MatrixXf tilde = X.transpose() * y;
  Eigen::VectorXf TrDinvXSX(k);
  Eigen::MatrixXf Dinv(p,k);
  if(TH){
    for(int i=0; i<k; i++){
      XSX.col(i) = XSX.col(i).array() * n(i);
    }
  }
  
  // Prior shape
  Eigen::MatrixXf Sb = vb*var_prior_df;
  Eigen::VectorXf Se = ve*var_prior_df;
  Eigen::VectorXf iNp = (n.array()+var_prior_df-1).inverse();
  
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  
  // Bending and convergence control
  Eigen::MatrixXf A = vb*1.0, GC(k,k);
  float bucketMean = 0.5*(bucketGCfrom+bucketGCto);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(A);
  float MinDVb, inflate = 0.0, Deflate = 1.0;
  Eigen::MatrixXf beta0(p,k), vb0(k,k);
  Eigen::VectorXf CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  float cnv = 10.0;
  int numit = 0;
  float logtol = log10(tol);
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Inner RGS
  std::vector<int> InnerRGSvec(k);
  for(int j=0; j<k; j++){InnerRGSvec[j]=j;}
  std::random_device rd2;
  std::mt19937 g2(rd2());
  int ri;
  
  // Non-Linear weights for marker effects
  bool NonLinear = NonLinearFactor!=0.0;
  Eigen::MatrixXf W(p,k);
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){  W(i,j) = 1.0; }}
  Eigen::VectorXf iVeWj = iVe*1.0;
  Eigen::VectorXf tmpW(p);
  float maxW, minW;
  
  // Objects for other variance structures
  float gs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(vb);
  Eigen::MatrixXf UDU(k,k);
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    vb0 = vb*1.0;
    ve0 = ve*1.0;
    h20 = h2*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    std::shuffle(InnerRGSvec.begin(), InnerRGSvec.end(), g2);
    
    for(int j=0; j<p; j++){
      
      J = RGSvec[j];
      b0 = b.row(J)*1.0;
      
      // System of equations - Traditional vs Stranden and Garrick 2009
      if(NoInversion){
        LHS = vb * (XX.row(J).transpose().array() * iVeWj.array()).matrix().asDiagonal(); 
        for(int i=0; i<k; i++){ LHS(i,i) += 1.0; }
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = (vb * (RHS.array() * iVeWj.array()).matrix()).array();
      }else{
        LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVeWj.array()).matrix();
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = RHS.array() * iVeWj.array();
      }
      
      // Inner GS
      if(InnerGS){
        b1 = b.row(J)*1.0;
        for(int i=0; i<k; i++){
          ri = InnerRGSvec[i];
          b1(ri) = (RHS(ri)-(LHS.col(ri).array()*b1.array()).sum()+LHS(ri,ri)*b1(ri))/LHS(ri,ri);}
      }else{
        b1 = LHS.llt().solve(RHS); 
      }
      
      // Update residuals
      b.row(J) = b1;
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Update marker weights
    if(NonLinear){
      W = b.cwiseAbs();
      for(int j=0; j<k; j++){
        maxW = W.col(j).maxCoeff(); minW = W.col(j).minCoeff();
        tmpW = NonLinearFactor * (W.col(j).array()-minW)/(maxW-minW) + (1.0-NonLinearFactor);
        tmpW = tmpW.array() + (1.0-tmpW.mean());
        W.col(j) = tmpW.array();
      }
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = (ve.array()+Se.array()) * iNp.array();
    h2 = 1 - ve.array()/vy.array();
    // Proportion-based prior
    if(weight_prior_h2>0){for(int i=0; i<k; i++){gs = ve(i)*(1-weight_prior_h2) + weight_prior_h2*veInit(i); ve(i) = gs*1.0;}}
    // Single variance
    if(OneVarE){tmp = ve.array().mean(); for(int i=0; i<k; i++) ve(i) = tmp*1.0;}
    iVe = ve.array().inverse();
    iVeWj = iVe*1.0;
    
    //Genetic variance
    
    // Get tilde-hat
    if(TH){
      for(int i=0; i<k; i++){
        Dinv.col(i) = (XSX.col(i).array()/ve(i) + iG(i,i)).inverse().array();
        TrDinvXSX(i)  = (XSX.col(i).transpose() * Dinv.col(i));
      }
      TildeHat = b.transpose()* Dinv.cwiseProduct(tilde);
    }else{
      TildeHat = b.transpose()*tilde;
    }
    
    // Estimate variances and covariance components
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ // Variances
          if(TH){
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrDinvXSX(i)+var_prior_df);
          }else{
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrXSX(i)+var_prior_df);
          }
        }else{ // Covariances
          if(TH){
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrDinvXSX(i)+TrDinvXSX(j)+var_prior_df);
          }else{
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrXSX(i)+TrXSX(j)+var_prior_df);
          }
        }}}
    
    
    if(weight_prior_h2>0){ // Proportion-based prior H2
      for(int i=0; i<k; i++){gs = vb(i,i)*(1-weight_prior_h2) + weight_prior_h2*vbInit(i); vb(i,i) = gs*1.0;}}
    if(weight_prior_gc>0){ // Proportion-based prior GC
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){
        if(i!=j){ GC(i,j) = (1.0-weight_prior_gc)*vb(i,j)/(sqrt(vb(i,i)*vb(j,j))) + gc_prior_df*weight_prior_gc;}else{GC(i,j) = 1.0;}}}
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){ if(i!=j){ vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}else{
        // Once calculation of GC without prior
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}}
      
      // Heterogeneous Compound Symmetry
      if(HCS){
        gs = 0.0;
        for(int i=0; i<k; i++){
          for(int j=0; j<k; j++){
            if(i>j){gs += GC(i,j);}}}
        gs = gs/((k*(k-1))/2);
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){ 
          if(i!=j){ GC(i,j) =  gs*1.0;}else{ GC(i,j) = 1.0; }}}
        // Extended Factor Analytics
      }else if(XFA){
        es.compute(GC);
        UDU = es.eigenvalues()[k] * es.eigenvectors().col(k) * es.eigenvectors().col(k).transpose();
        for(int i=1; i<NumXFA; i++) UDU += es.eigenvalues()[k-i] * es.eigenvectors().col(k-i) * es.eigenvectors().col(k-i).transpose();
        GC = UDU * 1.0; for(int i=0; i<k; i++){ GC(i,i)=1.0; };
      }
      
      // Monkeying with the correlations
      for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
          if(i!=j){
            // Zero'ing  Correlations
            if(MinCor<1.0){ if(GC(i,j)<MinCor){ GC(i,j) = 0.0; }}
            // Penalize Correlations
            if(PenCor>0.0){  GC(i,j) = tanh(PenCor*abs(GC(i,j)))*GC(i,j);} 
            // Round Down
            if(roundGCdownFrom<1.0){ if(GC(i,j)<roundGCdownFrom){ GC(i,j) = roundGCdownTo*1.0; }}
            // Round Up
            if(roundGCupFrom<1.0){ if(GC(i,j)>roundGCupFrom){ GC(i,j) = roundGCupTo*1.0; }}
            // Bucket round
            if(bucketGCfrom<1.0){ if(GC(i,j)>bucketGCfrom && GC(i,j)<bucketGCto  ){ GC(i,j) =  bucketMean*1.0; }}
            // Min H2
            if(uncorH2below>0.0){ if(h2(i)<uncorH2below || h2(j)<uncorH2below  ){ GC(i,j) = 0.0; }}
          }}}
      
      // RECONSTRUCT COVARIANCE HERE AND ONLY ONCE
      // Single variance of beta
      if(OneVarB){tmp = TildeHat.diagonal().mean(); vb = GC * tmp;  }else{
        // Regular covariance reconstruction
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){
          if(i!=j){ vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}
      
      // Bending
      if( !NoInversion || TH ){
        A = vb*1.0;
        for(int i=0; i<k; i++){ for(int j=0; j<k; j++){if(i!=j){A(i,j) *= Deflate;}}}
        if(A.llt().info()==Eigen::NumericalIssue && Deflate>DeflateMax){
          Deflate -= DeflateBy; if(verbose) Rcpp::Rcout <<  Deflate;
          for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*Deflate;}}}}
        EVDofA.compute(A);
        MinDVb = EVDofA.eigenvalues().minCoeff();
        if( MinDVb < 0.0 ){ 
          if(verbose) Rcpp::Rcout << ".";
          inflate = inflate - MinDVb*1.001;
          A.diagonal().array() += inflate;}
        iG = A.inverse();
      }
      
      // Compute convergence and print status
      
      //cnv = log10((beta0.array()-b.array()).square().sum());
      cnv = log10((beta0.array()-b.array()).square().colwise().sum().maxCoeff());
      CNV1(numit) = cnv; if(std::isnan(cnv)){ if(verbose){Rcpp::Rcout << "Numerical issue! Job aborted (it=" << numit << ")\n";} break;}
      CNV2(numit) = log10((h20.array()-h2.array()).square().sum());
      CNV3(numit) = log10((vb0.array()-vb.array()).square().sum());
      
      // Print
      ++numit;
      if( verbose && numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
      if( cnv<logtol ){ if(verbose){Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n";} break;}
      if( numit == maxit && verbose){ Rcpp::Rcout << "Model did not converge\n"; }
      
  }
  
  // Fitting the model
  Eigen::MatrixXf hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Resize convergence vectors
  Eigen::VectorXf CNV1b(numit),CNV2b(numit),CNV3b(numit);
  for(int i=0; i<numit; i++){
    CNV1b(i)=CNV1(i);
    CNV2b(i)=CNV2(i);
    CNV3b(i)=CNV3(i);
  }
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("vb")=vb,
                            Rcpp::Named("ve")=ve,
                            Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("bend")=Deflate,
                            Rcpp::Named("cnvB")=CNV1b,
                            Rcpp::Named("cnvH2")=CNV2b,
                            Rcpp::Named("cnvV")=CNV3b,
                            Rcpp::Named("b_Weights")=W,
                            Rcpp::Named("Its")=numit);
}
