#include <Rcpp.h>
using namespace Rcpp;

// Poisson-binomial mixture ---------------------------------------------------

// Calculate the density of a single draw from a
// Poisson-binomial mixture.
// This is an internal function only, and we assume
// the inputs are all correct.
// TODO - The normal approximation used by this function
// needs to be tested and refined.
// 'x' - The draw from the Poisson-binomial mixture
// 'size' - Sample size
// 'prob' - Probability
// 'use_log' - Whether to return the log density
// [[Rcpp::export]]
double dpoibin1(double x, double size, double prob, bool use_log) {
  double kThreshold = 50.0;
  double lambda = (1.0 - prob) * size;

  if (x > kThreshold) {
    double sd = sqrt((1.0 - pow(prob, 2.0)) * size);
    return R::dnorm(x, size, sd, use_log);
  }
  else {
    int limit = std::min(x, size);
    Range seq_binom = seq(0, limit);
    
    // seq doesn't allow 'inverse' sequences...
    Range seq_pois_needs_reversing = seq(x - limit, x); 
    // ...so we reverse the sequence order after creation
    IntegerVector seq_pois = rev(seq_pois_needs_reversing); 

    NumericVector logprob_binom = dbinom(seq_binom, size, prob, true);
    NumericVector logprob_pois = dpois(seq_pois, lambda, true);
    NumericVector logprob = logprob_binom + logprob_pois;
    
    if (use_log) {
      double max_val = max(logprob);
      return max_val + std::log(sum(exp(logprob - max_val)));
    }
    else {
      return sum(exp(logprob));
    }
  }
}


// Draw a single random variate from a Poisson-binomial
// mixture. This is an internal function only, and we assume
// the inputs are all correct.
// 'size' - Sample size
// 'prob' - Probability
// [[Rcpp::export]]
double rpoibin1(double size, double prob) {
  double val_binom = R::rbinom(size, prob);
  double val_pois = R::rpois((1.0 - prob) * size);
  return val_binom + val_pois;
}


// Poisson --------------------------------------------------------------------

// Draw a single random variate from a lower-truncated
// Poisson distribution.
// 'lambda' - the expected value (in the absence of truncation)
// 'lower' - the return value must greater than or
//   equal to 'lower'
// [[Rcpp::export]]
double rpoistr1(double lambda, double lower) {
  // first try rejection sampling...
  int max_attempt = 10;
  int n_attempt = 0;

  while (n_attempt++ < max_attempt) {
    double prop_value = R::rpois(lambda);
    if (prop_value >= lower)
      return prop_value;
  }

  //...and if that doesn't work, try sampling
  // based on cumulative distribution function 
  double p_lower = R::ppois(lower - 1.0, lambda, true, false);
  
  // If p_lower is too close to 1, then U
  // can be 1, which leads to a return
  // value of Inf. Note that 'p_lower_max' can't
  // be too close to 1 for this to work.
  constexpr double p_lower_max = 0.99999;

  if (p_lower > p_lower_max)
    return lower;

  double U = R::runif(p_lower, 1.0);
  
  double ans = R::qpois(U, lambda, true, false);

  // Due to rounding errors, 'ans' can
  // sometimes be lower than 'lower'. We
  // check and correct where necessary.
  if (ans < lower)
    return lower;

  // Also catch 'ans' that is effectively infinite,
  // which is easier than testing for actual infinity
  if (ans > 1e100)
    return lower;
  
  return ans;
}


// Draw 'n' values from lower-truncated Poisson distributions.
// This is an internal function only, and we assume the inputs
// are all correct.
// 'n' - the number of draws.
// 'lambda' - expected values (in the absence of truncation)
// 'lower' - each element of the return value must be
// greater than or equal to the corresponding value of 'lower'.
// [[Rcpp::export]]
NumericVector rpoistr(int n, NumericVector lambda, NumericVector lower) {
  NumericVector ans(n);

  for (int i=0; i < n; i++) {
    ans[i] = rpoistr1(lambda[i], lower[i]);
  }
  
  return ans;
}


// HAS_TESTS
// Calculate densitites for a truncated Poisson distribution.
// In places where 'x' is less than 'lower', return 0
// (if 'use_log' is false) or -Inf (if 'use_log' is true).
// 'x' - Draws from truncated Poisson distributions
// 'lambda' - Means for Poisson distributions
//    (in the absence of truncation). Same length as 'x'.
// 'lower' - All values of 'x' must be greater
//    than or equal to 'lower'. Same length as 'x'.
// 'use_log' - If 'true', return the log of the density.
// return value - Numeric vector the same length as 'x'.
// [[Rcpp::export]]
NumericVector dpoistr(NumericVector x,
		      NumericVector lambda,
		      NumericVector lower,
		      bool use_log) {
  int n = x.length();
  NumericVector ans(n);
  for (int i=0; i < n; i++) {
    if (x[i] < lower[i]) {
      ans[i] = use_log ? R_NegInf : 0;
    } else {
      double log_dens_untrunc = R::dpois(x[i], lambda[i], true);
      double log_const = R::ppois(lower[i] - 1, lambda[i], false, true);
      double val;
      if (is_infinite(NumericVector::create(log_dens_untrunc))[0]
	  || is_infinite(NumericVector::create(log_const))[0])
	val = R_NegInf;
      else
	val = log_dens_untrunc - log_const;
      ans[i] = use_log ? val : exp(val);
    }
  }
  return ans;
}

// pnchisq_approx -------------------------------------------------------------

// Fast approximation of cumulative distribution function
// for non-central chi-squared distribution, based on
// Sankaran 1959 (as reported in Wikipedia,
// https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution#Approximation_(including_for_quantiles)

// [[Rcpp::export]]
double pnchisq_approx(double q,
		      double df,
		      double ncp,
		      double lower,
		      bool use_log) {
  double A = df + ncp;
  double B = df + 2 * ncp;
  double C = df + 3 * ncp;
  double h = 1 - (2 * A * C) / (3 * B * B);
  double p = B / (A * A);
  double m = (h - 1) * (1 - 3 * h);
  double D = pow(q / A, h);
  double E = 1 + h * p * (h - 1 - 0.5 * (2 - h) * m * p);
  double F = h * sqrt(2 * p) * (1 + 0.5 * m * p);
  double v = (D - E) / F;
  return R::pnorm(v, 0, 1, lower, use_log);
}

// Switch to approximation when arguments fall outside
// bound set by 'threshold'

// [[Rcpp::export]]
double pnchisq_switch(double q,
		      double df,
		      double ncp,
		      double lower,
		      bool use_log) {
  constexpr double threshold = 20;
  if ((q > threshold) || (df > threshold) || (ncp > threshold))
    return pnchisq_approx(q, df, ncp, lower, use_log);
  else
    return R::pnchisq(q, df, ncp, lower, use_log);
}


// Skellam --------------------------------------------------------------------

// Calculate the cumulative distribution function
// for a quantile from a skellam distribution.
// When 'lower_tail' is TRUE, 'pskel' returns
// the probability that a skellam random variate
// with parameters 'mu1' and 'mu2' will be less
// than or equal to 'q'. When 'lower.tail' is
// FALSE, 'pskel' returns the probability
// that a skellam random variate with parameters
// 'mu1' and 'mu2' will be greater than 'q'.
// The calculations are modified from the function
// 'pskellam' in R package 'skellam', which
// in turn draw on 
// Johnson (1959) "On an Extension of the Connexion
// Between Poisson and Chi2 Distributions", Eq4.
// 'q' would normally be an integer, but
// the function can cope with non-integer values,
// so we allow these too.
// [[Rcpp::export]]
double pskel1(double q, double mu1, double mu2, bool lower_tail, bool log_p) {
  if (lower_tail) {
    if (q < 0)
      return pnchisq_switch(2.0 * mu2,
                           -2.0 * q,
                           2.0 * mu1,
                           true,
                           log_p);
      else
        return pnchisq_switch(2.0 * mu1,
                             2.0 * (q + 1.0),
                             2.0 * mu2,
                             false,
                             log_p);
  }
  else {
    if (q < 0)
      return pnchisq_switch(2.0 * mu2,
                           -2.0 * q,
                           2.0 * mu1,
                           false,
                           log_p);
      else
        return pnchisq_switch(2.0 * mu1,
                             2.0 * (q + 1.0),
                             2.0 * mu2,
                             true,
                             log_p);
  }
}


// Calculate the cumulative distribution function
// for quantities from a Skellam distribution.
// [[Rcpp::export]]
NumericVector pskel(NumericVector q,
		    NumericVector mu1,
		    NumericVector mu2,
		    bool lower_tail,
		    bool log_p) {
  int n = q.length();
  NumericVector ans(n);
  for (int i=0; i < n; i++) {
    ans[i] = pskel1(q[i], mu1[i], mu2[i], lower_tail, log_p);
  }
  return ans;
}


// Calculate value 'q' such that F(q) = p,
// where F is the cumulative distribution
// function for the Skellam distribution.
// This function is only ever used as a
// helper function for 'rskel1'. It always
// has 'lower.tail = TRUE' and 'log.p = FALSE'.
// Although 'q' is always an
// integer (in the mathematical sense)
// we represent it using a double.
// The function is modified from the function
// 'qskellam' in the 'skellam' package.
// [[Rcpp::export]]
double qskel1(double p, double mu1, double mu2) {
  // handle special cases where 'mu1' or 'mu2' is 0
  if (mu2 == 0.0) {
    return R::qpois(p,
                        mu1,
                        true,
                        false);
  }
  if (mu1 == 0.0) {
    return - R::qpois(p,
                         mu2,
                         false,
                         false);
  }
  
  // use the normal approximation, with
  // corrections for skew and kurtosis,
  // to get an initial value for 'q'
  double mu = mu1 - mu2;
  double sigma_sq = mu1 + mu2;
  double sigma = sqrt(sigma_sq);
  double z = R::qnorm(p,
                      0.0,
                      1.0,
                      true,
                      false);

  if (!R_finite(z))
    return z;

  double q_prop = mu + sigma * z;
  double skew = (pow(z, 2.0) - 1.0) * mu / sigma_sq / 6.0;
  double kurt = -(skew * mu - 2.0 * mu1 * mu2 * (pow(z, 2.0) - 3.0) / sigma_sq) * z / 12.0 / sigma;
  q_prop = round(q_prop + skew + kurt);

  // convert 'q_prop' to 'p_prop'
  double p_prop = pskel1(q_prop,
                   mu1,
                   mu2,
                   true,
                   false);
      
  // compare 'p_prop' to 'p' and if it is too high,
  // reduce 'q_prop', and hence 'p_prop' until
  // 'p_prop' is no longer too high
  bool too_high = p_prop > p;

  while (too_high) {
    q_prop--;
    p_prop = pskel1(q_prop,
                     mu1,
                     mu2,
                     true,
                     false);
    too_high = p_prop > p;
  }

  // compare 'p_prop' to 'p' and if it is too low,
  // increase 'q_prop', and hence 'p_prop' until
  // 'p_prop' is no longer too low
  bool too_low = p_prop < p;
  
  while (too_low) {
    q_prop++;
    p_prop = pskel1(q_prop,
                     mu1,
                     mu2,
                     true,
                     false);
    too_low = p_prop < p;
  }
    
  // value should now be just right, so return it
  return q_prop;
}



// Draw a single value from a left-truncated Skellam distribution.
// This is an internal function only, and we assume the inputs
// are all correct. We start by trying a simple rejection method,
// and if unsuccessful switch to the more reliable but slower
// inverse method. We guard against numerical problems
// by checking that 'lower' is not too far towards the right
// of the distribution.
// 'mu1', 'mu2' - mean paraketers
// 'lower' - the return value must
//    be greater than or equal to 'lower'.
// [[Rcpp::export]]
double rskeltr1(double mu1, double mu2, double lower) {

  // if 'lower' is higher than the 0.99999999th quantile,
  // then set then set the return value to 'lower'
  double upper_quant = qskel1(0.99999999, mu1, mu2);
  if (lower > upper_quant)
    return lower;
  
  // first try rejection sampling...
  int max_attempt = 10;
  int n_attempt = 0;

  while (n_attempt++ < max_attempt) {
    double x1 = R::rpois(mu1);
    double x2 = R::rpois(mu2);
    double prop_value = x1 - x2;

    if (prop_value >= lower)
      return prop_value;
  }

  // ...and if that doesn't work, try sampling
  // based on cumulative distribution function 
  double p_lower = pskel1(lower - 1.0, mu1, mu2, true, false);
  
  double U = R::runif(p_lower, 1.0);
  double ans = qskel1(U, mu1, mu2);
  
  // Due to rounding errors, 'ans' can
  // sometimes be lower than 'lower'. We
  // check and correct where necessary.
  if (ans < lower)
    return lower;
  
  return ans;
}

// Draw 'n' values from left-truncated Skellam distributions.
// This is an internal function only, and we assume the inputs
// are all correct.
// 'n' - a non-negative integer scalar. The number of draws.
// 'mu1', 'mu2' - mean parameters (in the absence of truncation)
// 'lower' - each element of
//    the return value must be greater than or equal to
//    the corresponding value of 'lower'.
// [[Rcpp::export]]
DoubleVector rskeltr(int n, DoubleVector mu1, DoubleVector mu2, DoubleVector lower) {
  DoubleVector ans(n);
  
  for (int i=0; i < n; i++) {
    ans[i] = rskeltr1(mu1[i], mu2[i], lower[i]);
  }
  
  return ans;
}


// Calculate the density of a draw from a
// Skellam distribution.
// This is an internal function only, and we assume
// the inputs are all correct.
// This function borrows the idea of using dpois when mu1 or mu2
// equals zero from the 'dskellam' function in package 'skellam'.
// It uses dnchisq rather than besselI (as skellam::dskellam
// does because dnchisq can cope with large values for its arguments.
// [[Rcpp::export]]
double dskel1(double x, double mu1, double mu2, bool use_log) {
  if (mu1 < 1e-6)
    return R::dpois(-x, mu2, use_log);
  if (mu2 < 1e-6)
    return R::dpois(x, mu1, use_log);
  double ans;
  if (x >= 0)
    ans = R::dnchisq(2 * mu1, 2 * (x + 1), 2 * mu2, use_log);
  else
    ans = R::dnchisq(2 * mu2, 2 * (1 - x), 2 * mu1, use_log);
  ans = use_log ? (log(2) + ans) : (2 * ans);
  return ans;
}
  

// Calculate the density of a draw from a
// left-truncated Skellam distribution.
// 'x' is drawn from a Skellam distribution,
// with the constraint that 'x' must be
// greater than or equal to 'lower'.
// The calculation of the density reflects the fact
// that 'rskeltr1' the quantile first.
// This is an internal function only, and we assume
// the inputs are all correct.
// [[Rcpp::export]]
double dskeltr1(double x,
		double mu1,
		double mu2,
		double lower,
		bool use_log) {
  if (x < lower)
    return (use_log ? R_NegInf : 0.0);

  // if 'lower' is higher than the 0.99999999th quantile,
  // then 'rskeltr1' would have set the return value to 'lower'
  double upper_quant = qskel1(0.99999999, mu1, mu2);
  if (lower > upper_quant)
    return (use_log ? 0 : 1);

  double log_dens_no_trunc = dskel1(x, mu1, mu2, true);
  double log_const = pskel1(lower - 1.0, mu1, mu2, false, true);
  double ans = log_dens_no_trunc - log_const;
  return (use_log ? ans : exp(ans));
}

// Calculate densities from a lower-truncated
// Skellam distribution.
// [[Rcpp::export]]
NumericVector dskeltr(NumericVector x,
		      NumericVector mu1,
		      NumericVector mu2,
		      NumericVector lower,
		      bool use_log) {
  int n = x.length();
  NumericVector ans(n);
  for (int i = 0; i < n; i++) {
    ans[i] = dskeltr1(x[i], mu1[i], mu2[i], lower[i], use_log);
  }
  return ans;
}




