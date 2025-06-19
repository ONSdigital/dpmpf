
#include <Rcpp.h>
using namespace Rcpp;

// update_counts --------------------------------------------------------------

// Copy values from 'stock_end', 'counts_bth', ..., 'counts_em2'
// to 'counts_stock', 'counts_births', ..., 'counts_emigration2'

// [[Rcpp::export]]
void update_counts_noreg_inner(Environment &self, int i_interval) {
  int n_particle = self["n_particle"];
  NumericVector counts_stock(self["counts_stock"]);
  NumericVector counts_births(self["counts_births"]);
  NumericVector counts_deaths(self["counts_deaths"]);
  NumericVector counts_immigration1(self["counts_immigration1"]);
  NumericVector counts_emigration1(self["counts_emigration1"]);
  NumericVector counts_immigration2(self["counts_immigration2"]);
  NumericVector counts_emigration2(self["counts_emigration2"]);
  NumericVector stock_end(self["stock_end"]);
  NumericVector counts_bth(self["counts_bth"]);
  NumericVector counts_dth(self["counts_dth"]);
  NumericVector counts_im1(self["counts_im1"]);
  NumericVector counts_em1(self["counts_em1"]);
  NumericVector counts_im2(self["counts_im2"]);
  NumericVector counts_em2(self["counts_em2"]);
  int offset_stock = (i_interval + 1) * n_particle;
  int offset_event = i_interval * n_particle;
  for (int i = 0; i < n_particle; i++) {
    int i_stock = offset_stock + i;
    int i_event = offset_event + i;
    counts_stock[i_stock] = stock_end[i];
    counts_births[i_event] = counts_bth[i];
    counts_deaths[i_event] = counts_dth[i];
    counts_immigration1[i_event] = counts_im1[i];
    counts_emigration1[i_event] = counts_em1[i];
    counts_immigration2[i_event] = counts_im2[i];
    counts_emigration2[i_event] = counts_em2[i];
  }
  LogicalVector is_parallelogram_first(self["is_parallelogram_first"]);
  if (is_parallelogram_first[i_interval]) {
    NumericVector stock_end_second(self["stock_end_second"]);
    NumericVector counts_bth_second(self["counts_bth_second"]);
    NumericVector counts_dth_second(self["counts_dth_second"]);
    NumericVector counts_im1_second(self["counts_im1_second"]);
    NumericVector counts_em1_second(self["counts_em1_second"]);
    int offset_stock = (i_interval + 2) * n_particle;
    int offset_event = (i_interval + 1) * n_particle;
    for (int i = 0; i < n_particle; i++) {
      int i_stock = offset_stock + i;
      int i_event = offset_event + i;
      counts_stock[i_stock] = stock_end_second[i];
      counts_births[i_event] = counts_bth_second[i];
      counts_deaths[i_event] = counts_dth_second[i];
      counts_immigration1[i_event] = counts_im1_second[i];
      counts_emigration1[i_event] = counts_em1_second[i];
    }
  }
}

// [[Rcpp::export]]
void update_counts_withreg_inner(Environment &self, int i_interval) {
  int n_particle = self["n_particle"];
  int n_region = self["n_region"];
  NumericVector counts_stock(self["counts_stock"]);
  NumericVector counts_births(self["counts_births"]);
  NumericVector counts_deaths(self["counts_deaths"]);
  NumericVector counts_internal_in(self["counts_internal_in"]);
  NumericVector counts_internal_out(self["counts_internal_out"]);
  NumericVector counts_immigration1(self["counts_immigration1"]);
  NumericVector counts_emigration1(self["counts_emigration1"]);
  NumericVector counts_immigration2(self["counts_immigration2"]);
  NumericVector counts_emigration2(self["counts_emigration2"]);
  NumericVector stock_end(self["stock_end"]);
  NumericVector counts_bth(self["counts_bth"]);
  NumericVector counts_dth(self["counts_dth"]);
  NumericVector counts_in(self["counts_in"]);
  NumericVector counts_out(self["counts_out"]);
  NumericVector counts_im1(self["counts_im1"]);
  NumericVector counts_em1(self["counts_em1"]);
  NumericVector counts_im2(self["counts_im2"]);
  NumericVector counts_em2(self["counts_em2"]);
  int length_slice = n_particle * n_region;
  int offset_stock = (i_interval + 1) * length_slice;
  int offset_event = i_interval * length_slice;
  for (int i = 0; i < length_slice; i++) {
    int i_stock = offset_stock + i;
    int i_event = offset_event + i;
    counts_stock[i_stock] = stock_end[i];
    counts_births[i_event] = counts_bth[i];
    counts_deaths[i_event] = counts_dth[i];
    counts_internal_in[i_event] = counts_in[i];
    counts_internal_out[i_event] = counts_out[i];
    counts_immigration1[i_event] = counts_im1[i];
    counts_emigration1[i_event] = counts_em1[i];
    counts_immigration2[i_event] = counts_im2[i];
    counts_emigration2[i_event] = counts_em2[i];
  }
}
