#include <Rcpp.h>
using namespace Rcpp;


// Log of sum of two exponentials
//
// Log of sum of exponentials, calculated
// in a way that avoids overflow/underflow
// Assume 'x' and 'y' have equal length.
// [[Rcpp::export]]
NumericVector log_sum_exp_2(NumericVector x, NumericVector y) {
  int n = x.length();
  NumericVector ans(n);
  for (int i = 0; i < n; i++) {
    double xi = x[i];
    double yi = y[i];
    bool both_neginf = (is_infinite(NumericVector::create(xi))[0]
			&& is_infinite(NumericVector::create(yi))[0]
			&& (xi < 0)
			&& (yi < 0));
    if (both_neginf)
      ans[i] = R_NegInf;
    else if (xi > yi)
      ans[i] = xi + log1p(exp(yi - xi));
    else
      ans[i] = yi + log1p(exp(xi - yi));
  }
  return ans;	
}
