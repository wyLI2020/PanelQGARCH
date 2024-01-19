#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double fc_ifCover(double observation, double q) {
  double ifCover = 0.0;
  if(observation < q) {
    ifCover = 1.0;
  } else {
    ifCover = 0.0;
  }
  return ifCover;
}

//[[Rcpp::export]]
arma::vec fc_ifCovers_vec(arma::vec observations_vec, arma::vec q_vec) {
  int N = observations_vec.n_elem;
  arma::mat ifCovers_mat(N, 1); ifCovers_mat.fill(0.0);
  for(int i = 0; i < N; i++) {
    ifCovers_mat(i, 0) = fc_ifCover(observations_vec(i), q_vec(i));
  }
  arma::vec ifCovers_vec = ifCovers_mat.col(0);
  return ifCovers_vec;
}

// double fc_ECR(arma::rowvec observations_vec, double q) {
//   int t1 = observations_vec.n_elem;
//   arma::mat ifCovers_mat(t1, 1); ifCovers_mat.fill(0.0);
//   for(int t = 0; t < t1; t++) {
//     ifCovers_mat(t, 0) = fc_ifCover(observations_vec(t), q);
//   }
//   arma::vec ifCovers = ifCovers_mat.col(0);
//   double ECR = mean(ifCovers);
//   return ECR;
// }
// 
// //[[Rcpp::export]]
// double fc_average_ECR(arma::mat Y_N_t1, arma::vec q_N) {
//   int N = Y_N_t1.n_rows;
//   arma::mat ECRs_mat(N, 1); ECRs_mat.fill(0.0);
//   for(int i = 0; i < N; i++) {
//     ECRs_mat(i, 0) = fc_ECR(Y_N_t1.row(i), q_N(i));
//   }
//   arma::vec ECRs = ECRs_mat.col(0);
//   double average_ECR = mean(ECRs);
//   return average_ECR;
// }

