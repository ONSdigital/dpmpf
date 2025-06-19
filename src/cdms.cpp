#include <Rcpp.h>
#include "cdm.h"
using namespace Rcpp;


// 'Cdms' classes are lists of cohort data models. 


class CdmsNoreg {
public:
  CdmsNoreg(List models) : models(models) {};
  
  NumericVector calc_loglik(NumericVector& counts_true,
			    int i_interval,
			    double obs_zero) {
    int n_particle = counts_true.length();
    NumericVector ans (n_particle);
    int n_model = models.length();
    
    for (int i=0; i < n_model; i++) {
      CdmNoregBase& model = as<CdmNoregBase&>(models[i]);

      ans += model.calc_loglik(counts_true, i_interval, obs_zero);
    }
    
    return ans;
  }
  
  NumericVector draw_counts_true(int i_interval, int n_particle) {
    NumericVector counts_true = rep(NA_REAL, n_particle);
    int n_model = models.length();

    for (int i=0; i < n_model; i++) {
      CdmNoregBase& model = as<CdmNoregBase&>(models[i]);

      model.fill_counts_true(counts_true,
                     i_interval);

      bool is_empty = R_IsNA(counts_true[0]);

      if (!is_empty)
        break;
    }

    return counts_true;
  }
  
  NumericVector calc_logimp(NumericVector& counts_true,
			    int i_interval) {
    int n_particle = counts_true.length();
    NumericVector logimp = rep(NA_REAL, n_particle);
    int n_model = models.length();

    for (int i=0; i < n_model; i++) {
      CdmNoregBase& model = as<CdmNoregBase&>(models[i]);

      model.fill_logimp(logimp, counts_true, i_interval);

      bool is_empty = R_IsNA(logimp[0]);
      
      if (!is_empty)
        break;
    }

    return logimp;
  }

  List models;
};

RCPP_EXPOSED_CLASS(CdmsNoreg)
  RCPP_MODULE(mod_cdmsnoreg) {
    class_<CdmsNoreg>("CdmsNoreg")
    .constructor<List>("cohort-level data models, no region")
    .method( "calc_loglik", &CdmsNoreg::calc_loglik )
    .method( "draw_counts_true", &CdmsNoreg::draw_counts_true )
    .method( "calc_logimp", &CdmsNoreg::calc_logimp )
    ;
  }


class CdmsWithreg {
public:
  CdmsWithreg(List models) : models(models) {};

  NumericVector calc_loglik(NumericMatrix& counts_true,
			    int i_interval,
			    double obs_zero) {
    int n_particle = counts_true.rows();
    NumericVector ans (n_particle);
    int n_model = models.length();
    
    for (int i=0; i < n_model; i++) {
      CdmWithregBase& model = as<CdmWithregBase&>(models[i]);
      
      ans += model.calc_loglik(counts_true, i_interval, obs_zero);
    }
    
    return ans;
  }

  NumericMatrix draw_counts_true(int i_interval,
				 int n_particle,
				 int n_region) {

    NumericMatrix counts_true (n_particle, n_region);
    std::fill(counts_true.begin(), counts_true.end(), NA_REAL);

    int n_model = models.length();

    for (int i=0; i < n_model; i++) {
      CdmWithregBase& model = as<CdmWithregBase&>(models[i]);
      
      model.fill_counts_true(counts_true,
			     i_interval);
      
      NumericMatrix::Row first_particle = counts_true(0,_);
      
      bool has_empty = false;
      for(auto value : first_particle){
	if(R_IsNA(value)) {
	  has_empty = true;
	  break;
	}
      }
      
      if (!has_empty)
	break;
    }

    return counts_true;
  }

  NumericMatrix calc_logimp(NumericMatrix& counts_true,
                            int i_interval) {
    int n_particle = counts_true.rows();
    int n_region = counts_true.cols();
    NumericMatrix logimp (n_particle, n_region);
    std::fill(logimp.begin(), logimp.end(), NA_REAL);
    
    int n_model = models.length();
    
    for (int i=0; i < n_model; i++) {
      CdmWithregBase& model = as<CdmWithregBase&>(models[i]);
      
      model.fill_logimp(logimp, counts_true, i_interval);
      
      NumericMatrix::Row first_particle = logimp(0,_);
      
      bool has_empty = false;
      for(auto value : first_particle){
        if(R_IsNA(value)) {
          has_empty = true;
          break;
        }
      }
      
      if (!has_empty)
        break;
    }
    
    return logimp;
  }
  
  List models;
};

RCPP_EXPOSED_CLASS(CdmsWithreg)
  RCPP_MODULE(mod_cdmswithreg) {
    class_<CdmsWithreg>("CdmsWithreg")
    .constructor<List>("cohort-level data models, with region")
    .method( "calc_loglik", &CdmsWithreg::calc_loglik )
    .method( "draw_counts_true", &CdmsWithreg::draw_counts_true )
    .method( "calc_logimp", &CdmsWithreg::calc_logimp )
    ;
  }
