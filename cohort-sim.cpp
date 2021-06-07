// [[Rcpp::depends(BH)]]

#include <boost/multi_array.hpp>
#include <Rcpp.h>
#include <fstream>
#include <iostream>

using namespace Rcpp;
using namespace boost;

// [[Rcpp::export]]
NumericVector cohort_sim(const List s_fp) {
  const List s_ss = s_fp["ss"];
  const int n_year(s_ss["PROJ_YEARS"]), n_sex(2), n_age(66), n_agr(9), n_cd4(7);
  const int n_ts(10);
  const double dt(0.1);

  // Map single age indices (0=15..65=80) to age groups
  const int a2g[n_age] = {
    0, 0, 1, 1, 1, 2, 2, 2, 2, 2, // 15-16 = 0, 17-19 = 1, 20-24 = 2
    3, 3, 3, 3, 3, 4, 4, 4, 4, 4, // 25-29 = 3, 30-34 = 4
    5, 5, 5, 5, 5, 6, 6, 6, 6, 6, // 35-39 = 5, 40-44 = 6
    7, 7, 7, 7, 7, 8, 8, 8, 8, 8, // 45-49 = 7, 50-54 = 8
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, // 55-64 = 8
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, // 65-74 = 8
    8, 8, 8, 8, 8, 8};            // 75-80 = 8

  multi_array_ref<double, 3> surv(REAL(s_fp["Sx"]), extents[n_year][n_sex][n_age]);
  multi_array_ref<double, 3> inci(REAL(s_fp["new_inf"]), extents[n_year][n_sex][n_age]);
  multi_array_ref<double, 3> in_dist(REAL(s_fp["cd4_initdist"]), extents[n_sex][n_agr][n_cd4]);
  multi_array_ref<double, 3> in_prog(REAL(s_fp["cd4_prog"]), extents[n_sex][n_agr][n_cd4-1]);
  multi_array_ref<double, 3> in_mort(REAL(s_fp["cd4_mort"]), extents[n_sex][n_agr][n_cd4]);
  
  multi_array<double, 3> dist(extents[n_sex][n_age][n_cd4]);
  multi_array<double, 3> prog(extents[n_sex][n_age][n_cd4-1]);
  multi_array<double, 3> mort(extents[n_sex][n_age][n_cd4]);
  
  double hivpop[n_year * n_sex * n_age * n_cd4];
  multi_array_ref<double, 4> Y(hivpop, extents[n_year][n_sex][n_age][n_cd4]);
  double influx[n_cd4], efflux[n_cd4];

  memset(hivpop, 0, n_year * n_sex * n_age * n_cd4 * sizeof(double));

  for (int si(0); si < n_sex; ++si) {
    for (int ai(0); ai < n_age; ++ai) {
      for (int hi(0); hi < n_cd4; ++hi) dist[si][ai][hi] = in_dist[si][a2g[ai]][hi];
      for (int hi(0); hi < n_cd4; ++hi) mort[si][ai][hi] = in_mort[si][a2g[ai]][hi] * dt; // 1.0 - exp(-in_mort[si][a2g[ai]][hi] * dt);
      for (int hi(0); hi < n_cd4-1; ++hi) prog[si][ai][hi] = in_prog[si][a2g[ai]][hi] * dt; // 1.0 - exp(-in_prog[si][a2g[ai]][hi] * dt);
    }
  }
  
  // Year 1 simulation
  for (int si(0); si < n_sex; ++si) {
    for (int ai(0); ai < n_age; ++ai) {
      for (int hi(0); hi < n_cd4; ++hi) {
        Y[0][si][ai][hi] = inci[0][si][ai] * dist[si][ai][hi];
      }
    }
  }

  // Subsequent year simulation
  for (int yi(1); yi < n_year; ++yi) {

    // Aging and non-HIV mortality
    for (int si(0); si < n_sex; ++si) {
      for (int hi(0); hi < n_cd4; ++hi) {
        for (int ai(1); ai < n_age; ++ai) Y[yi][si][ai][hi] = Y[yi-1][si][ai-1][hi] * surv[yi][si][ai];
        Y[yi][si][n_age-1][hi] += Y[yi-1][si][n_age-1][hi] * surv[yi][si][n_age-1];
      }
    }

    // HIV disease progression and mortality
    for (int ti(0); ti < n_ts; ++ti) {
      for (int si(0); si < n_sex; ++si) {
        for (int ai(0); ai < n_age; ++ai) {
          efflux[n_cd4-1] = Y[yi][si][ai][n_cd4-1] * mort[si][ai][n_cd4-1];
          for (int hi(0); hi < n_cd4-1; ++hi) efflux[hi] = Y[yi][si][ai][hi] * (prog[si][ai][hi] + mort[si][ai][hi]);
          for (int hi(1); hi < n_cd4; ++hi) influx[hi] = Y[yi][si][ai][hi-1] * prog[si][ai][hi-1];
          for (int hi(0); hi < n_cd4; ++hi) Y[yi][si][ai][hi] += influx[hi] - efflux[hi];
        }
      }
    }
    
    // Distribute new HIV infections
    for (int si(0); si < n_sex; ++si) {
      for (int ai(0); ai < n_age; ++ai) {
        for (int hi(0); hi < n_cd4; ++hi) Y[yi][si][ai][hi] += inci[yi][si][ai] * dist[si][ai][hi];
      }
    }
  }
  
  NumericVector rval(n_year * n_sex * n_age * n_cd4);
  for (int i(0); i < rval.length(); ++i) rval[i] = hivpop[i];
  return rval; // in R, convert to array via array(rval, c(n_cd4, n_age, n_sex, n_year))
}
