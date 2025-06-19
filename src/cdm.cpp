#include "cdm.h"
#include "distributions.h"

CdmNoregPoibin::CdmNoregPoibin(NumericVector counts_data, double prob) : counts_data(counts_data), prob(prob) {}

NumericVector CdmNoregPoibin::calc_loglik(NumericVector& counts_true,
					  int i_interval,
					  double obs_zero) {
  int n_particle = counts_true.length();
  double val_data  = counts_data[i_interval];
  NumericVector ans = rep(0.0, n_particle);

  if (!R_IsNA(val_data)) {
    for (int i_particle = 0; i_particle < n_particle; i_particle++) {
      double val_true = counts_true[i_particle];
      bool is_zero = (val_true < 1);
      if (is_zero)
	ans[i_particle] = R::dpois(val_data,
				   obs_zero,
				   true);
      else
	ans[i_particle] = dpoibin1(val_data,
				   val_true,
				   prob,
				   true);
    }
  }
  return ans;
}

void CdmNoregPoibin::fill_counts_true(NumericVector& counts_true,
				      int i_interval) {
  double val_data = counts_data[i_interval];
  
  if (R_IsNA(val_data))
    return;
    
  int n_particle = counts_true.length();

  for (int i_particle = 0; i_particle < n_particle; i_particle++) {
    counts_true[i_particle] = rpoibin1(val_data,
				       prob);
  }
}
  
void CdmNoregPoibin::fill_logimp(NumericVector& logimp,
				 NumericVector& counts_true,
				 int i_interval) {
  double val_data = counts_data[i_interval];
  
  if (R_IsNA(val_data))
    return;
    
  int n_particle = logimp.length();

  for (int i = 0; i < n_particle; i++) {
    double val_true = counts_true[i];
    logimp[i] = dpoibin1(val_true,
			 val_data,
			 prob,
			 true);
  }
}

  

RCPP_MODULE(mod_cdmnoregpoibin) {
  class_<CdmNoregPoibin>("CdmNoregPoibin")
    .constructor<NumericVector, double>("Cohort data model based on Poisson-binomial mixture, with no regions")
    .method( "calc_loglik", &CdmNoregPoibin::calc_loglik )
    .method( "fill_counts_true", &CdmNoregPoibin::fill_counts_true )
    .method( "fill_logimp", &CdmNoregPoibin::fill_logimp )
  ;
}


CdmWithregPoibin::CdmWithregPoibin(NumericMatrix counts_data, double prob) : counts_data(counts_data), prob(prob) {}

NumericVector CdmWithregPoibin::calc_loglik(NumericMatrix& counts_true,
					    int i_interval,
					    double obs_zero) {
  int n_particle = counts_true.rows();
  int n_region = counts_true.cols();
  NumericVector ans (n_particle);

  for (int i_region = 0; i_region <  n_region; i_region++) {
    double val_data = counts_data(i_region, i_interval);

    if (!R_IsNA(val_data)) {
      for (int i_particle = 0; i_particle < n_particle; i_particle++) {
        double val_true = counts_true(i_particle, i_region);
	bool is_zero = (val_true < 1);
	if (is_zero)
	  ans[i_particle] += R::dpois(val_data,
				      obs_zero,
				      true);
	else
	  ans[i_particle] += dpoibin1(val_data,
				      val_true,
				      prob,
				      true);
      }
    }
  }

  return ans;
}

void CdmWithregPoibin::fill_counts_true(NumericMatrix& counts_true,
					int i_interval) {
  int n_particle = counts_true.rows();
  int n_region = counts_true.cols();
  
  // some columns of 'counts_true' can be filled
  // while others are empty, because previous datasets
  // had observations for some regions but not others
  for (int i_region = 0; i_region < n_region; i_region++) {
    double val_true = counts_true(0, i_region);
    bool region_has_been_filled = !R_IsNA(val_true);
    
    if (region_has_been_filled)
      continue;
    
    double val_data = counts_data(i_region, i_interval);
    bool have_data_for_region = !R_IsNA(val_data);
    
    if (!have_data_for_region)
        continue;
    
    for (int i_particle = 0; i_particle < n_particle; i_particle++) {
      counts_true(i_particle, i_region) = rpoibin1(val_data, prob);
    }
  }
}
  
void CdmWithregPoibin::fill_logimp(NumericMatrix& logimp,
				   NumericMatrix& counts_true,
				   int i_interval) {
  int n_particle = counts_true.rows();
  int n_region = counts_true.cols();

  for (int i_region = 0; i_region < n_region; i_region++) {
    double val_logimp = logimp(0, i_region);
    bool region_has_been_filled = !R_IsNA(val_logimp);

    if (region_has_been_filled)
      continue;

    double val_data = counts_data(i_region, i_interval);
    bool have_data_for_region = !R_IsNA(val_data);
    if (!have_data_for_region)
      continue;

    for (int i_particle = 0; i_particle < n_particle; i_particle++) {
      double val_true = counts_true(i_particle, i_region);
      logimp(i_particle, i_region) = dpoibin1(val_true,
					      val_data,
					      prob,
					      true);
    }
  }
}


RCPP_MODULE(mod_cdmwithregpoibin) {
  class_<CdmWithregPoibin>("CdmWithregPoibin")
  .constructor<NumericMatrix, double>("Cohort data model based on Poisson-binomial mixture, with regions")
  .field( "counts_data", &CdmWithregPoibin::counts_data )
  .method( "calc_loglik", &CdmWithregPoibin::calc_loglik )
  .method( "fill_counts_true", &CdmWithregPoibin::fill_counts_true )
  .method( "fill_logimp", &CdmWithregPoibin::fill_logimp )
  ;
}



CdmNoregNbinom::CdmNoregNbinom(NumericVector counts_data,
                               NumericVector ratio,
                               NumericVector disp) : counts_data(counts_data), ratio(ratio), disp(disp) {}

NumericVector CdmNoregNbinom::calc_loglik(NumericVector& counts_true,
					  int i_interval,
					  double obs_zero) {
  constexpr double threshold_zero = 0.000001;

  int n_particle = counts_true.length();
  double val_data  = counts_data[i_interval];
  double val_ratio = ratio[i_interval];
  double val_disp = disp[i_interval];
  bool disp_is_zero = val_disp < threshold_zero;

  NumericVector ans = rep(0.0, n_particle);

  if (!R_IsNA(val_data)) {
    for (int i_particle = 0; i_particle < n_particle; i_particle++) {
      double val_true = counts_true[i_particle];
      double mu = val_ratio * val_true;
      bool val_is_zero = (val_true < 1.0);
      if (val_is_zero) {
	ans[i_particle] = R::dpois(val_data, obs_zero, true);
      } else if (disp_is_zero) {
	ans[i_particle] = R::dpois(val_data, mu, true);
      } else {
	double size = 1.0 / val_disp;
	ans[i_particle] = R::dnbinom_mu(val_data, size, mu, true);
      }
    }
  }

  return ans;
}


void CdmNoregNbinom::fill_counts_true(NumericVector& counts_true,
				      int i_interval) {
  constexpr double threshold_zero = 0.000001;
  double val_data = counts_data[i_interval];
  
  if (R_IsNA(val_data))
    return;
  
  // can't derive a value, so leave blank
  double val_ratio = ratio[i_interval];
  if (val_ratio < threshold_zero)
    return;

  double mu = val_data / val_ratio;
  double val_disp = disp[i_interval];
  bool disp_is_zero = val_disp < threshold_zero;

  int n_particle = counts_true.length();
  for (int i_particle = 0; i_particle < n_particle; i_particle++) {
    if (disp_is_zero) {
      counts_true[i_particle] = R::rpois(mu);
    } else {
      double size = 1.0 / val_disp;
      double prob = size / (size + mu);
      counts_true[i_particle] = R::rnbinom(size, prob); // rnbinom_mu not working
    }
  }
}


void CdmNoregNbinom::fill_logimp(NumericVector& logimp,
                                 NumericVector& counts_true,
                                 int i_interval) {
  constexpr double threshold_zero = 0.000001;
  double val_data = counts_data[i_interval];
  
  if (R_IsNA(val_data))
    return;
  
  double val_ratio = ratio[i_interval];
  if (val_ratio < threshold_zero)
    return;
  
  double mu = val_data / val_ratio;
  double val_disp =  disp[i_interval];
  bool disp_is_zero = val_disp < threshold_zero;
  
  int n_particle = logimp.length();
  for (int i = 0; i < n_particle; i++) {
    double val_true = counts_true[i];
    if (disp_is_zero) {
      logimp[i] = R::dpois(val_true, mu, true);
    } else {
      double val_size = 1.0 / val_disp;
      logimp[i] = R::dnbinom_mu(val_true, val_size, mu, true);
    }
  }
}

RCPP_MODULE(mod_cdmnoregnbinom) {
  class_<CdmNoregNbinom>("CdmNoregNbinom")
  .constructor<NumericVector, NumericVector, NumericVector>("Cohort data model based on negative binomial distribution, with no regions")
  .field( "counts_data", &CdmNoregNbinom::counts_data )
  .method( "calc_loglik", &CdmNoregNbinom::calc_loglik )
  .method( "fill_counts_true", &CdmNoregNbinom::fill_counts_true )
  .method( "fill_logimp", &CdmNoregNbinom::fill_logimp )
  ;
}



CdmWithregNbinom::CdmWithregNbinom(NumericMatrix counts_data, NumericMatrix ratio, NumericMatrix disp) :
  counts_data(counts_data), ratio(ratio), disp(disp) {}

NumericVector CdmWithregNbinom::calc_loglik(NumericMatrix& counts_true,
					    int i_interval,
					    double obs_zero) {
  constexpr double threshold_zero = 0.000001;

  int n_particle = counts_true.rows();
  int n_region = counts_true.cols();
  NumericVector ans (n_particle);

  for (int i_region = 0; i_region <  n_region; i_region++) {
    double val_data = counts_data(i_region, i_interval);
    if (!R_IsNA(val_data)) {
      double val_ratio = ratio(i_region, i_interval);
      double val_disp = disp(i_region, i_interval);
      bool disp_is_zero = val_disp < threshold_zero;

      double val_size = 1.0 / disp(i_region, i_interval);
      for (int i_particle = 0; i_particle < n_particle; i_particle++) {
        double val_true = counts_true(i_particle, i_region);
        double mu = val_ratio * val_true;
	bool val_is_zero = (val_true < 1.0);
	if (val_is_zero) {
	  ans[i_particle] += R::dpois(val_data, obs_zero, true);
	} else if (disp_is_zero) {
	  ans[i_particle] += R::dpois(val_data, mu, true);
	} else {
	  ans[i_particle] += R::dnbinom_mu(val_data, val_size, mu, true);
	}
      }
    }
  }
  
  return ans;
}

void CdmWithregNbinom::fill_counts_true(NumericMatrix& counts_true,
					int i_interval) {
  constexpr double threshold_zero = 0.000001;
  int n_particle = counts_true.rows();
  int n_region = counts_true.cols();
  
  // some columns of 'counts_true' can be filled
  // while others are empty, because previous datasets
  // had observations for some regions but not others
  for (int i_region = 0; i_region < n_region; i_region++) {
    double val_true = counts_true(0, i_region);

    bool region_has_been_filled = !R_IsNA(val_true);

    if (region_has_been_filled)
      continue;

    double val_data = counts_data(i_region, i_interval);
    bool have_data_for_region = !R_IsNA(val_data);
    
    if (!have_data_for_region)
      continue;

    // can't derive a value, so leave blank
    double val_ratio = ratio(i_region, i_interval);
    if (val_ratio < threshold_zero)
    return;

    double val_disp = disp(i_region, i_interval);
    bool disp_is_zero = val_disp < threshold_zero;
    
    for (int i_particle = 0; i_particle < n_particle; i_particle++) {
      double mu = val_data / val_ratio;
      if (disp_is_zero) {
	counts_true(i_particle, i_region) = R::rpois(mu);
      } else {
	double size = 1.0 / val_disp;
	double prob = size / (size + mu);
	counts_true(i_particle, i_region) = R::rnbinom(size, prob); // rnbinom_mu not working
      }
    }
  }
}

void CdmWithregNbinom::fill_logimp(NumericMatrix& logimp,
                                   NumericMatrix& counts_true,
                                   int i_interval) {
  constexpr double threshold_zero = 0.000001;
  int n_particle = counts_true.rows();
  int n_region = counts_true.cols();
  
  for (int i_region = 0; i_region < n_region; i_region++) {
    double val_logimp = logimp(0, i_region);
    bool region_has_been_filled = !R_IsNA(val_logimp);
    
    if (region_has_been_filled)
      continue;
    
    double val_data = counts_data(i_region, i_interval);
    bool have_data_for_region = !R_IsNA(val_data);
    if (!have_data_for_region)
      continue;
    
    double val_ratio = ratio(i_region, i_interval);
    if (val_ratio < threshold_zero)
      return;
    
    double mu = val_data / val_ratio;
    double val_disp = disp(i_region, i_interval);
    bool disp_is_zero = val_disp < threshold_zero;
    
    for (int i_particle = 0; i_particle < n_particle; i_particle++) {
      double val_true = counts_true(i_particle, i_region);
      if (disp_is_zero) {
	logimp(i_particle, i_region) = R::dpois(val_true, mu, true);
      } else {
	double size = 1.0 / val_disp;
	logimp(i_particle, i_region) = R::dnbinom_mu(val_true, size, mu, true);
      }
    }
  }
}

RCPP_MODULE(mod_cdmwithregnbinom) {
  class_<CdmWithregNbinom>("CdmWithregNbinom")
  .constructor<NumericMatrix, NumericMatrix, NumericMatrix>("Cohort data model baysed on negative binomial distribution, with regions")
  .field( "counts_data", &CdmWithregNbinom::counts_data )
  .method( "calc_loglik", &CdmWithregNbinom::calc_loglik )
  .method( "fill_counts_true", &CdmWithregNbinom::fill_counts_true )
  .method( "fill_logimp", &CdmWithregNbinom::fill_logimp )
  ;
}


CdmNoregNorm::CdmNoregNorm(NumericVector counts_data,
                           NumericVector ratio,
                           NumericVector sd) : counts_data(counts_data), ratio(ratio), sd(sd) {}

NumericVector CdmNoregNorm::calc_loglik(NumericVector& counts_true,
					  int i_interval,
					  double obs_zero) {

  int n_particle = counts_true.length();
  double val_data  = counts_data[i_interval];
  double val_ratio = ratio[i_interval];
  double val_sd = sd[i_interval];

  NumericVector ans = rep(0.0, n_particle);

  if (!R_IsNA(val_data)) {
    for (int i_particle = 0; i_particle < n_particle; i_particle++) {
      double val_true = counts_true[i_particle];
      double mu = val_ratio * val_true;
      bool val_is_zero = (val_true < 1.0);
      if (val_is_zero) {
	ans[i_particle] = R::dpois(val_data, obs_zero, true);
      } else {
	double upper = R::pnorm(val_data + 0.5, mu, val_sd, true, true);
	double lower = R::pnorm(val_data - 0.5, mu, val_sd, true, true);
	ans[i_particle] = upper + log1p(-exp(lower - upper));
      }
    }
  }

  return ans;
}


void CdmNoregNorm::fill_counts_true(NumericVector& counts_true,
				      int i_interval) {
  constexpr double threshold_zero = 0.000001;
  double val_data = counts_data[i_interval];
  
  if (R_IsNA(val_data))
    return;
  
  // can't derive a value, so leave blank
  double val_ratio = ratio[i_interval];
  if (val_ratio < threshold_zero)
    return;

  double mu = val_data / val_ratio;
  double val_sd = sd[i_interval];

  int n_particle = counts_true.length();
  for (int i_particle = 0; i_particle < n_particle; i_particle++) {
    counts_true[i_particle] = static_cast<int>(R::rnorm(mu, val_sd) + 0.5);
  }
}


void CdmNoregNorm::fill_logimp(NumericVector& logimp,
                                 NumericVector& counts_true,
                                 int i_interval) {
  constexpr double threshold_zero = 0.000001;
  double val_data = counts_data[i_interval];
  
  if (R_IsNA(val_data))
    return;
  
  double val_ratio = ratio[i_interval];
  if (val_ratio < threshold_zero)
    return;
  
  double mu = val_data / val_ratio;
  double val_sd =  sd[i_interval];
  
  int n_particle = logimp.length();
  for (int i = 0; i < n_particle; i++) {
    double val_true = counts_true[i];
    double upper = R::pnorm(val_true + 0.5, mu, val_sd, true, true);
    double lower = R::pnorm(val_true - 0.5, mu, val_sd, true, true);
    logimp[i] = upper + log1p(-exp(lower - upper));
  }
}



RCPP_MODULE(mod_cdmnoregnorm) {
  class_<CdmNoregNorm>("CdmNoregNorm")
  .constructor<NumericVector, NumericVector, NumericVector>("Cohort data model based on normal distribution, with no regions")
  .field( "counts_data", &CdmNoregNorm::counts_data )
  .method( "calc_loglik", &CdmNoregNorm::calc_loglik )
  .method( "fill_counts_true", &CdmNoregNorm::fill_counts_true )
  .method( "fill_logimp", &CdmNoregNorm::fill_logimp )
  ;
}
