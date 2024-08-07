#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computePL(NumericMatrix D, int homDim, NumericVector scaleSeq,int k=1){
  int n_rows = 0; // number of rows with the correct dimension
  for(int i=0;i<D.nrow();++i){
    if(D(i,0) == homDim){
      ++n_rows; 
    }
  }
  
  int L = scaleSeq.size();
  if (n_rows == 0) return NumericVector(L);
  
  NumericVector x(n_rows),y(n_rows);
  int n=0;
  for(int i=0;i<D.nrow();++i){
    if(D(i,0) == homDim){
      x[n] = D(i,1);
      y[n] = D(i,2);
      ++n;
    }
  }
  
  NumericVector lambda(L);
  for (int i=0;i<L;++i){
    NumericVector Lambda = pmax(pmin(scaleSeq[i] - x, y - scaleSeq[i]),0);
    lambda[i] = Lambda.sort(true)[k-1];
  }
  return lambda; 
}

