#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computePS(NumericMatrix D, int homDim, NumericVector scaleSeq,int p=1){
  int n_rows = 0; // number of rows with the correct dimension
  for(int i=0;i<D.nrow();++i){
    if((D(i,0) == homDim)&&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
      ++n_rows; 
    }
  }
  
  if (n_rows == 0) return NumericVector(scaleSeq.size()-1);
  
  NumericVector x(n_rows),y(n_rows);
  int n=0;
  for(int i=0;i<D.nrow();++i){
    if((D(i,0) == homDim)&&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
      x[n] = D(i,1);
      y[n] = D(i,2);
      ++n;
    }
  }
  
  int l = scaleSeq.size()-1; 
  NumericVector pp = pow(y-x,p);
  NumericVector w = pp/sum(pp);
  
  NumericVector phi(l);
  for (int k=0;k<l;++k){
    NumericVector alpha1 = pmax(scaleSeq[k],x);
    NumericVector alpha2 = pmax(scaleSeq[k],(x+y)/2);
    
    NumericVector beta1 = pmin(scaleSeq[k+1],(x+y)/2);
    NumericVector beta2 = pmin(scaleSeq[k+1],y);
    
    NumericVector b1 = pmax(0,beta1-alpha1)*((beta1+alpha1)/2-x);
    NumericVector b2 = pmax(0,beta2-alpha2)*(y-(beta2+alpha2)/2);
    phi[k] = sum(w*(b1+b2))/(scaleSeq[k+1]-scaleSeq[k]);
  }
  return phi;  
}

