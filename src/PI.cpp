#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

NumericVector seq_C(double a, double b, double by){
  int n = round((b-a)/by)+1;
  NumericVector out(n);
  for(int i=0; i<n; ++i){
    out[i] = a + i*by;
  }
  return out;
}

NumericVector outer_C(NumericVector x, NumericVector y) {
  int n = x.length();
  int m = y.length();
  NumericVector result(n*m);
  for ( int i = 0; i < m; ++i ) {
    for ( int j = 0; j < n; ++j ) {
      result[j+i*n] = x[j] * y[i];
    }
  }
  return result;
}

// PS for H0
NumericVector PSurfaceH0(NumericVector point,
                         NumericVector y_lower,
                         NumericVector y_upper,
                         double sigma,
                         double maxP) {
  double y = point[1];
  NumericVector out2 = pnorm(y_upper,y,sigma) - pnorm(y_lower,y,sigma);
  double wgt = y/maxP * (y < maxP) + (y >= maxP);
  return wgt*out2;
}

// PS for Hk, k>0
NumericVector PSurfaceHk(NumericVector point,
                         NumericVector y_lower,
                         NumericVector y_upper,
                         NumericVector x_lower,
                         NumericVector x_upper,
                         double sigma,
                         double maxP) {
  double x = point[0];
  double y = point[1];
  NumericVector out1 = pnorm(x_upper,x,sigma) - pnorm(x_lower,x,sigma);
  NumericVector out2 = pnorm(y_upper,y,sigma) - pnorm(y_lower,y,sigma);
  double wgt = y/maxP * (y < maxP) + (y >= maxP);
  return wgt*outer_C(out1,out2); 
}

// [[Rcpp::export]]
NumericVector computePI(NumericMatrix D,int homDim,
                        NumericVector xSeq, NumericVector ySeq,
                        double sigma){
// D - N by 3 matrix (columns contain dimension, birth and persistence values respectively)
int n_rows = 0; // number of rows with the correct dimension
  for(int i=0; i<D.nrow(); ++i) {
    if((D(i,0) == homDim)&&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
      ++n_rows; 
    }
  }
  

int resP = ySeq.size()-1; 
int resB = xSeq.size()-1;
if (n_rows == 0){
  if (resB>0){
    return NumericVector(resB*resP);
  } 
  else{
    return NumericVector(resP);
  }
}

NumericMatrix D_(n_rows,2);
int j=0;
for(int i=0; i<D.nrow(); ++i) {
  if((D(i,0) == homDim)&&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
    D_(j,0) = D(i,1);
    D_(j,1) = D(i,2);
    ++j;
  }
}

double minP = ySeq[0];
double maxP = ySeq[ySeq.size()-1];
double dy = (maxP-minP)/resP;
NumericVector y_lower = seq_C(minP,maxP-dy,dy);
NumericVector y_upper = y_lower + dy;
int Nsize;
double sumB = sum(abs(diff(D_(_,0))));
if ((homDim==0)&&(sumB==0)) {
  Nsize = resP;
  }else{
  Nsize = resB*resP;
}
NumericMatrix Psurf_mat(Nsize,n_rows);
  
if (Nsize==resP) {
  for(int i=0; i<n_rows; ++i){
    Psurf_mat(_,i) = PSurfaceH0(D_(i,_),y_lower,y_upper,sigma,maxP);
  }

} else{
  double minB = xSeq[0];
  double maxB = xSeq[xSeq.size()-1];
  double dx = (maxB-minB)/resB;
  NumericVector x_lower = seq_C(minB,maxB-dx,dx);
  NumericVector x_upper = x_lower + dx;
  for(int i=0; i<n_rows; ++i){
    Psurf_mat(_,i) = PSurfaceHk(D_(i,_),y_lower,y_upper,x_lower,x_upper,sigma,maxP);
  }
}
NumericVector out = rowSums(Psurf_mat);
return out;
}
