#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Self-weight function based on He & Yi (2021)
//'
//'
//'
//' @param Y Vector. Data.
//' @param C Double. Quantile of Y or |Y|.
//'
//' @export
NumericVector weight_function_HeYi2021_cpp(NumericVector Y,double C){
  int N=Y.length(),t,i;
  NumericVector w(N);

  for(t=0;t<N;t++){
    for(i=1;i<=t;i++){
      w[t]+=exp(-pow(log(i),2))*((fabs(Y[t-i])<C)? 1:(fabs(Y[t-i])/C));
    }

    for(i=t+1;i<N;i++){
      w[t]+=exp(-pow(log(i),2))*1;
    }
    w[t] = pow(w[t],-3);
  }

  return w;
}

//[[Rcpp::export]]
arma::mat fc_weighs_NT(arma::mat Y_N_T, arma::vec c_N) {
  // matrix of weights, i.e. [w_{it}] for i=1,...,N and t=1,...,T
  int N = Y_N_T.n_rows;
  int T = Y_N_T.n_cols;
  arma::mat t_Y_T_N = Y_N_T.t();
  arma::mat t_w_T_N(T, N); t_w_T_N.fill(0.0);
  for(int i = 0; i < N; i++) {
    arma::vec y_i_T = t_Y_T_N.col(i);
    double c_i = c_N(i);
    NumericVector y_i_T_nv(y_i_T.begin(), y_i_T.end());
    NumericVector w_i_T_nv = weight_function_HeYi2021_cpp(y_i_T_nv, c_i);
    arma::colvec w_i_T(w_i_T_nv.begin(), w_i_T_nv.size(), false);
    t_w_T_N.col(i) = w_i_T;
  }
  arma::mat w_N_T = t_w_T_N.t();
  return w_N_T;
}

arma::vec fc_phitilde_i_tau(double tau, arma::vec y_i_T, arma::mat Factors_r_T, arma::vec w_i_T) {
  // \widetilde{\bm{\phi}}_{i}(\tau)
  Environment OPT = Environment::namespace_env("PanelQGARCHwithFactor");
  Function f = OPT["fr_phitilde_i_tau"];
  // Function f("fr_phitilde_i_tau");
  NumericVector phitilde_i_nv = f(tau, y_i_T, Factors_r_T, w_i_T);
  arma::vec phitilde_i(phitilde_i_nv.begin(), phitilde_i_nv.size(), false);
  return phitilde_i;
}

//[[Rcpp::export]]
arma::mat fc_phitilde_N_tau(double tau, arma::mat Y_N_T, arma::mat Factors_r_T, arma::mat w_N_T) {
  // the N \times 3+r matrix of initial estimates with the i-th row being \widetilde{\bm{\phi}}_{i}(\tau)
  int N = Y_N_T.n_rows;
  int r = Factors_r_T.n_rows;
  arma::mat phitilde_3plusr_N(3+r, N); phitilde_3plusr_N.fill(0.0);
  arma::mat t_Y_T_N = Y_N_T.t();
  arma::mat t_w_T_N = w_N_T.t();
  for(int i = 0; i < N; i++) {
    arma::vec y_i_T = t_Y_T_N.col(i);
    arma::vec w_i_T = t_w_T_N.col(i);
    arma::vec phitilde_i = fc_phitilde_i_tau(tau, y_i_T, Factors_r_T, w_i_T);
    phitilde_3plusr_N.col(i) = phitilde_i;
  }
  arma::mat phitilde_N_3plusr = phitilde_3plusr_N.t(); 
  return phitilde_N_3plusr;
}

// //[[Rcpp::export]]
// arma::vec fc_q_y_t0plus1(arma::mat Y_N_t0, arma::mat phihat_N_3plusr) {
//   // (Q_\tau(y_{1, t_0 + 1}|\mathcal{F}_{N-1}), ..., Q_\tau(y_{N, t_0 + 1}|\mathcal{F}_{N-1}))
//   Function f("q_y");
//   int N = Y_N_t0.n_rows;
//   int t0 = Y_N_t0.n_cols;
//   arma::mat Y_t0_N = Y_N_t0.t();
//   arma::mat phihat_3plusr_N = phihat_N_3plusr.t();
//   arma::mat q_y_t0plus1_mat(N, 1); q_y_t0plus1_mat.fill(0.0);
//   for(int i = 1; i < (N+1); i++) {
//     arma::vec y_i = Y_t0_N.col(i-1);
//     arma::vec phihat_i = phihat_3plusr_N.col(i-1);
//     NumericVector y_i_nv(y_i.begin(), y_i.end());
//     NumericVector phihat_i_nv(phihat_i.begin(), phihat_i.end());
//     NumericVector q_y_i_2TOt0plus1_nv = f(phihat_i_nv, y_i_nv);
//     arma::vec q_y_i_2TOt0plus1(q_y_i_2TOt0plus1_nv.begin(), q_y_i_2TOt0plus1_nv.size(), false);
//     double q_y_i_t0plus1 = q_y_i_2TOt0plus1(t0-1);
//     q_y_t0plus1_mat(i-1, 0) = q_y_i_t0plus1;
//   }
//   arma::vec q_y_t0plus1 = q_y_t0plus1_mat.col(0);
//   return q_y_t0plus1;
// }

