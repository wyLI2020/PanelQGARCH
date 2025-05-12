#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec fc_partition_phi_l_N(int l, arma::vec g_vec, arma::vec changepoints_vec, arma::mat rank_mat) {
  // obtain the partition of (\bm\phi_{l1}, ..., \bm\phi_{lN})' from the change points and the ranks, when l = 1,2,3, ..., 3+r
  Environment OPT = Environment::namespace_env("PanelQGARCHwithFactor");
  Function f = OPT["fr_partition_phi_l_N"];
  // Function f("fr_partition_phi_l_N");
  NumericVector partition_phi_l_N_nv = f(l, g_vec, changepoints_vec, rank_mat);
  arma::vec partition_phi_l_N(partition_phi_l_N_nv.begin(), partition_phi_l_N_nv.size(), false);
  return partition_phi_l_N;
}

// [[Rcpp::export]]
arma::mat fc_partition_phi_N(arma::vec g_vec, arma::vec changepoints_vec, arma::mat rank_mat) {
  // obtain the partition of (\bm\phi_1, ..., \bm\phi_N)' from the change points and the ranks
  int N = rank_mat.n_rows; 
  int rplus3 = g_vec.n_elem; 
  arma::mat partition_phi_N_3plusr(N,rplus3); partition_phi_N_3plusr.fill(0.0); 
  for (int l = 1; l < (rplus3+1); l++) {
    partition_phi_N_3plusr.col(l-1) = fc_partition_phi_l_N(l, g_vec, changepoints_vec, rank_mat);
  }
  return(partition_phi_N_3plusr); 
}

