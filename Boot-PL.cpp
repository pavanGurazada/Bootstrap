#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

/**
 * @title Bootstrap from a power law distribution
 * @author Pavan Gurazada
 * @licence MIT
 * @summary Generate bootstrap samples and fit a power law distribution
 *
*/


double alphaPLFit(const arma::vec& x) {
  Environment igraph("package:igraph");
  Function fit_power_law = igraph["fit_power_law"];
  
  List p = fit_power_law(Named("x", x));
  
  return p["alpha"];
}


double xMinPLFit(const arma::vec& x) {
  Environment igraph("package:igraph");
  Function fit_power_law = igraph["fit_power_law"];
  
  List p = fit_power_law(Named("x", x));
  
  return p["xmin"];
}

// [[Rcpp::export]]
arma::mat bootstraps(const arma::vec& x, int times) {
  arma::mat bootMatrix;
  
  for (int t = 0; t < times; ++t) {
    arma::vec resample = RcppArmadillo::sample(x, x.n_elem, true);
    bootMatrix.insert_cols(t, resample);
  }
  
  return bootMatrix;
}

// [[Rcpp::export]]
arma::vec bootAlphaPLFit(const arma::vec& x, int times) {
  arma::mat bootMatrix = bootstraps(x, times);
  arma::vec alphas(times);
  
  for (int i = 0; i < times; ++i) {
    double alpha = alphaPLFit(bootMatrix.col(i));
    alphas[i] = alpha;
  }
  
  return alphas;
}

// [[Rcpp::export]]
arma::vec bootxMinPLFit(const arma::vec& x, int times) {
  arma::mat bootMatrix = bootstraps(x, times);
  arma::vec xMins(times);
  
  for (int i = 0; i < times; ++i) {
    double xMin = xMinPLFit(bootMatrix.col(i));
    xMins[i] = xMin;
  }
  
  return xMins;
}

/***R

*/