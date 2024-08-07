#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeNL(NumericMatrix D, int homDim, NumericVector scaleSeq){
  int n_rows = 0; // number of rows with the correct dimension
  int scaleLen = scaleSeq.size()-1;
  for(int i=0;i<D.nrow();++i){
    if((D(i,0) == homDim)&&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
      ++n_rows; 
    }
  }
  
  if (n_rows == 0) return NumericVector(scaleLen);
  
  NumericVector x(n_rows),y(n_rows);
  int n=0;
  for(int i=0;i<D.nrow();++i){
    if((D(i,0) == homDim)&&(Rcpp::traits::is_finite<REALSXP>(D(i,2)))){
      x[n] = D(i,1);
      y[n] = D(i,2);
      ++n;
    }
  }
 
  NumericVector lL = (y - x)/sum(y-x); //compute weights

  NumericVector nl(scaleLen);
  NumericVector b(n);
  for (int k=0;k<scaleLen;++k){
    b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x);
    nl[k] = sum(lL*pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]);
  }
  return nl; 
}

