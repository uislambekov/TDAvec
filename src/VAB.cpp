#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeVAB(NumericMatrix D, int homDim, NumericVector scaleSeq){
  int n_rows = 0; // number of rows with the correct dimension
  for(int i=0;i<D.nrow();++i){
    if(D(i,0) == homDim){
      ++n_rows; 
    }
  }
  
  if (n_rows == 0) return NumericVector(scaleSeq.size()-1);
  
  NumericVector x(n_rows),y(n_rows);
  int n=0;
  for(int i=0;i<D.nrow();++i){
    if( D(i,0) == homDim){
      x[n] = D(i,1);
      y[n] = D(i,2);
      ++n;
    }
  }
  
  int l = scaleSeq.size()-1; 
  NumericVector vab(l);
  for (int k=0;k<l;++k){
    NumericVector b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x);
    vab[k] = sum(pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]);
  }
  return vab; 
}

