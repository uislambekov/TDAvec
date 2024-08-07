#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

IntegerVector which_C(LogicalVector vec) {
  int n = sum(vec);
  IntegerVector outC(n);
  int k =0;
  for(int i=0; i<vec.size(); ++i) {
    if(vec[i]) {
      outC[k]=i;
      ++k;
    }
  }
  return outC;
}

// [[Rcpp::export]]
NumericVector computeVPB(NumericMatrix D, int homDim,
                                NumericVector xSeq, NumericVector ySeq,
                                double tau=0.3) {
  IntegerVector which_C(LogicalVector vec);
  int n_rows = 0; // number of rows with the correct dimension
  for(int i=0; i<D.nrow(); ++i) {
    if(D(i,0) == homDim) {
      n_rows = n_rows+1;
    }
  }

  NumericVector x(n_rows), y(n_rows);
  int j=0;
  for(int i=0; i<D.nrow(); ++i) {
    if( D(i,0) == homDim) {
      x[j] = D(i,1);
      y[j] = D(i,2);
      ++j;
    }
  }

  NumericVector lambda = tau*y;
  NumericVector dy = diff(ySeq);
  int m = dy.size();
  double sumB = sum(abs(diff(x)));
  if ((homDim==0)&&(sumB==0)) {
    NumericVector vpb(m);
    for(int i=0; i<m; ++i) {
      double c = ySeq[i];
      double d = ySeq[i+1];
      vpb[i] = 0;
      for(int j=0; j<y.size(); ++j) {
        if( (y[j] > c - lambda[j]) && (y[j] < d+lambda[j])) {
          double y_cd = y[j];
          double lambda_cd = lambda[j];
          double yMin = std::max(c, y_cd - lambda_cd);
          double yMax = std::min(d, y_cd + lambda_cd);
          vpb[i] += 0.5*(yMax*yMax-yMin*yMin)/dy[i] ;
        }
      }
    }
    return vpb;
  }
  else{
    NumericVector dx = diff(xSeq);
    int n = dx.size();
    NumericVector vpb(n*m);
    for(int i=1; i<=m; ++i) {
      double c = ySeq[i-1];
      double d =  ySeq[i+1-1];
      IntegerVector yInd =which_C( (y>c-lambda) & (y<d+lambda));
      if (yInd.size()>0){
        NumericVector y_cd = y[yInd];
        NumericVector x_cd = x[yInd];
        NumericVector lambda_cd = lambda[yInd];
        NumericVector yMin_cd = pmax(c,y_cd-lambda_cd);
        NumericVector yMax_cd = pmin(d,y_cd+lambda_cd);
        for (int j=1; j<=n; ++j) {
          double a = xSeq[j-1];
          double b = xSeq[j+1-1];
          IntegerVector xInd = which_C((x_cd>a-lambda_cd) & (x_cd<b+lambda_cd));
          if (xInd.size()>0) {
            NumericVector x_abcd = x_cd[xInd];
            NumericVector lambda_abcd = lambda_cd[xInd];
            NumericVector xMin = pmax(a,x_abcd-lambda_abcd);
            NumericVector xMax = pmin(b,x_abcd+lambda_abcd);
            NumericVector yMin = yMin_cd[xInd];
            NumericVector yMax = yMax_cd[xInd];
            vpb[j+(i-1)*n-1] = 0.5*sum((xMax-xMin)*(yMax-yMin)*(xMin+xMax+yMin+yMax))/(dx[j-1]*dy[i-1]);
          }
          else vpb[j+(i-1)*n-1] = 0;
        }
      }
    }
    return vpb;
  }
}


