#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

/**
 * @title Simulating a first order autoregressive process
 * @author Pavan Gurazada
 * @licence MIT
 * @summary https://www.r-bloggers.com/another-nice-rcpp-example/
 *
*/

// [[Rcpp::export]]
arma::mat simulate(arma::mat a, arma::mat e) {
  arma::mat coef = a;
  arma::mat errors = e;
  
  int m = errors.n_rows;
  int n = errors.n_cols;
  arma::mat sim_data(m, n);
  
  sim_data.row(0) = arma::zeros(1, n);
  
  for (int r = 1; r < m; ++r) {
    sim_data.row(r) = sim_data.row(r-1) * arma::trans(coef) + errors.row(r);
  }
  
  return sim_data;
  
}

/***R
a <- matrix(c(0.5, 0.1, 0.1, 0.5), nrow = 2)
e <- matrix(rnorm(1e4), ncol = 2)

data <- simulate(a, e)

*/