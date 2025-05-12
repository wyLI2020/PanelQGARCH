#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec fc_rank(arma::vec vec_1) {
  // R function "rank"
  Environment OPT = Environment::namespace_env("base");
  Function f = OPT["rank"];
  // Function f("rank");
  NumericVector rank_vec_nv = f(vec_1);
  arma::vec rank_vec(rank_vec_nv.begin(), rank_vec_nv.size(), false);
  return rank_vec;
}

// [[Rcpp::export]]
arma::mat fc_rank_mat(arma::mat phi_N_3plusr) {
  // the matrix of ranks for \widetilde\bm\phi
  int N = phi_N_3plusr.n_rows; 
  int rplus3 = phi_N_3plusr.n_cols; 
  arma::mat rank_mat(N, rplus3); rank_mat.fill(0.0); 
  for (int l = 0; l < rplus3; l++) {
    rank_mat.col(l) = fc_rank(phi_N_3plusr.col(l));
  }
  return rank_mat;
}

arma::vec fc_phi_l_N_homo(int l, arma::vec theta_vec, arma::vec g_vec, arma::vec changepoints_vec, arma::mat rank_mat) {
  // transform the parameter vector (\omega_{(1)}, ..., \omega_{(g_1)})', (\alpha_{(1)}, ..., \alpha_{(g_2)})', (\beta_{(1)}, ..., \beta_{(g_3)})', or (\lambda_{l(1)}, ..., \lambda_{l(g_{3+r})})' 
  // to the N \times 1 parameter vector (\omega_{1}, ..., \omega_{N})', (\alpha_{1}, ..., \alpha_{N})', (\beta_{1}, ..., \beta_{N})', or (\lambda_{l1}, ..., \lambda_{lN})' under the homogeneity structure
  Environment OPT = Environment::namespace_env("PanelQGARCHwithFactor");
  Function f = OPT["fr_phi_l_N_homo"];
  // Function f("fr_phi_l_N_homo");
  NumericVector phi_l_N_nv = f(l, theta_vec, g_vec, changepoints_vec, rank_mat);
  arma::vec phi_l_N(phi_l_N_nv.begin(), phi_l_N_nv.size(), false);
  return phi_l_N;
}

// [[Rcpp::export]]
arma::mat fc_phi_N_homo(arma::vec theta_vec, arma::vec g_vec, arma::vec changepoints_vec, arma::mat rank_mat) {
  // transform the parameter vector \bm\theta = (\omega_{(1)}, ..., \omega_{(g_1)}, \alpha_{(1)}, ..., \alpha_{(g_2)}, \beta_{(1)}, ..., \beta_{(g_3)}, \lambda_{1(1)}, ..., \lambda_{1(g_4)}, ..., \lambda_{r(1)}, ..., \lambda_{r(g_{3+r})})' 
  // to the N \times (3+r) parameter matrix (\bm{\phi}_1, ..., \bm{\phi}_N)' under the homogeneity structure
  int N = rank_mat.n_rows; 
  int rplus3 = g_vec.n_elem; 
  arma::mat phi_N_3plusr(N,rplus3); phi_N_3plusr.fill(0.0); 
  for (int l = 1; l < (rplus3+1); l++) {
    phi_N_3plusr.col(l-1) = fc_phi_l_N_homo(l, theta_vec, g_vec, changepoints_vec, rank_mat);
  }
  return phi_N_3plusr; 
}

// [[Rcpp::export]]
arma::vec fc_calculate_segment_means(arma::vec vector_ascend, arma::vec points) {
  // calculate the mean of each segment, where "points" are the split points of "vector_ascend" with the fist point being 0 and the last point being the length of "vector_ascend"
  arma::vec means_vec(points.n_elem - 1); 
  for (int i = 0; i < points.n_elem - 1; i++) {
    int start = points(i); //        0,   a1, a2, ...
    int end = points(i+1) - 1; // a1-1, a2-1,     ..., N-1
    means_vec(i) = mean(vector_ascend.subvec(start, end)); 
  }
  return means_vec;
}

arma::vec fc_theta_l_initial_IN_optim(int l, arma::mat phiini_N_3plusr, arma::vec g_vec, arma::vec changepoints_vec) {
  // initial value in optim() for estimating \widehat\bm\omega, \widehat\bm\alpha, \widehat\bm\beta, or \widehat\bm\lambda_l
  Environment OPT = Environment::namespace_env("PanelQGARCHwithFactor");
  Function f = OPT["fr_theta_l_initial_IN_optim"];
  // Function f("fr_theta_l_initial_IN_optim");
  NumericVector theta_l_ini_nv = f(l, phiini_N_3plusr, g_vec, changepoints_vec);
  arma::vec theta_l_ini(theta_l_ini_nv.begin(), theta_l_ini_nv.size(), false);
  return theta_l_ini;
}

// [[Rcpp::export]]
arma::vec fc_theta_initial_IN_optim(arma::mat phiini_N_3plusr, arma::vec g_vec, arma::vec changepoints_vec) {
  // initial value in optim() for estimating \widehat\bm\theta
  int rplus3 = g_vec.n_elem; 
  arma::vec theta_ini(sum(g_vec)); 
  int g_1 = g_vec(0);
  arma::vec theta_1_ini = fc_theta_l_initial_IN_optim(1, phiini_N_3plusr, g_vec, changepoints_vec);
  theta_ini.subvec(0, g_1-1) = theta_1_ini;
  for (int l = 2; l < (rplus3+1); l++) {
    int g_l = g_vec(l-1);
    arma::vec theta_l_ini = fc_theta_l_initial_IN_optim(l, phiini_N_3plusr, g_vec, changepoints_vec);
    theta_ini.subvec(sum(g_vec.subvec(0,l-2)), sum(g_vec.subvec(0,l-2))+g_l-1) = theta_l_ini;
  }
  return theta_ini;
}

// [[Rcpp::export]]
double fc_Delta_ijh(int i, int j, int h, arma::vec vec_ascend) {
  // \Delta_{i,j}(h) for 1 \leq i \leq N
  arma::vec vec_i_h = vec_ascend.subvec(i-1, h-1);
  arma::vec vec_hplus1_j = vec_ascend.subvec(h, j-1);
  double average_i_h = mean(vec_i_h);
  double average_hplus1_j = mean(vec_hplus1_j);
  double Delta_ijh = sqrt(1.0 * (j-h) * (h-i+1) / (j-i+1)) * std::abs(average_hplus1_j - average_i_h);
  return Delta_ijh;
}

int fc_hhat_maxDelta_ij(int i, int j, arma::vec vec_ascend) {
  // \widehat{h}_{i,j} = max_{i \leq h < j} \Delta_{i,j}(h)
  arma::vec Deltas_ijh(j-i); Deltas_ijh.fill(0.0);
  for(int h = i; h < j; h++) {
    Deltas_ijh(h-i) = fc_Delta_ijh(i, j, h, vec_ascend);
  }
  int hhat = index_max(Deltas_ijh) + i;
  return hhat;
}

// [[Rcpp::export]]
arma::Col<int> fc_changepoints(double delta, int i, int j, arma::vec vec_ascend) {
  // obtain the change points among (i,j) using the binary segmentation algorithm, including the point j
  arma::Mat<int> hhats_mat(1,1); hhats_mat(0,0) = j;
  if(j-i+1 < 2) { return hhats_mat.col(0); }
  int hhat_ij = fc_hhat_maxDelta_ij(i, j, vec_ascend);
  double Delta_hhat_ij = fc_Delta_ijh(i, j, hhat_ij, vec_ascend);
  if(Delta_hhat_ij < delta || Delta_hhat_ij == delta) { return hhats_mat.col(0); }
  arma::Mat<int> hhat_ij_mat(1,1); hhat_ij_mat(0,0) = hhat_ij; 
  hhats_mat = join_vert(hhats_mat, hhat_ij_mat);
  hhats_mat = join_vert(hhats_mat, fc_changepoints(delta, i, hhat_ij, vec_ascend));
  hhats_mat = join_vert(hhats_mat, fc_changepoints(delta, hhat_ij+1, j, vec_ascend));
  return sort(unique(hhats_mat.col(0)));
}

// [[Rcpp::export]]
arma::Col<int> fc_g_changepoints_vec(arma::mat phi_N_3plusr, arma::vec delta_123r) {
  // the vector of the numbers of change points and change points for \widetilde\bm\phi
  int N = phi_N_3plusr.n_rows; 
  int rplus3 = phi_N_3plusr.n_cols; 
  arma::Col<int> g_vec(rplus3);
  double delta_1 = delta_123r(0); // \widehat{\delta}_1
  arma::vec phi_1_tilde_N = phi_N_3plusr.col(0); // when l=1, (\widetilde\omega_{1}, ..., \widetilde\omega_{N})'
  arma::vec phi_1_tilde_N_ascend = sort(phi_1_tilde_N); // when l=1, (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  arma::Col<int> changepoints_1 = fc_changepoints(delta_1, 1, N, phi_1_tilde_N_ascend); // the vector of change points when l=1
  g_vec(0) = changepoints_1.n_elem;
  arma::Col<int> changepoints_vec = changepoints_1; 
  for (int l = 2; l < (rplus3+1); l++) {
    double delta_l = delta_123r(l-1); 
    arma::vec phi_l_tilde_N = phi_N_3plusr.col(l-1); 
    arma::vec phi_l_tilde_N_ascend = sort(phi_l_tilde_N); 
    arma::Col<int> changepoints_l = fc_changepoints(delta_l, 1, N, phi_l_tilde_N_ascend); 
    g_vec(l-1) = changepoints_l.n_elem;
    changepoints_vec = join_cols(changepoints_vec, changepoints_l);
  }
  arma::Col<int> g_changepoints_vec = join_cols(g_vec, changepoints_vec);
  return g_changepoints_vec;
}

double fc_loss_QR(double tau, arma::vec phi_i, arma::vec y_i_T, arma::mat Factors_r_T, arma::vec w_i_T) {
  // the loss function of the initial estimator for individual i
  Environment OPT = Environment::namespace_env("PanelQGARCHwithFactor");
  Function f = OPT["Loss_QR"];
  // Function f("Loss_QR");
  NumericVector loss_i_nv = f(phi_i, y_i_T, Factors_r_T, w_i_T, tau);
  arma::vec loss_i(loss_i_nv.begin(), loss_i_nv.size(), false);
  return loss_i(0);
}

// [[Rcpp::export]]
double fc_loss_panelQR(double tau, arma::mat phi_N_3plusr, arma::mat Y_N_T, arma::mat Factors_r_T, arma::mat w_N_T) {
  // the loss function of the final estimator 
  int N = Y_N_T.n_rows;
  double loss_PanelQR = 0.0;
  arma::mat t_phi_3plusr_N = phi_N_3plusr.t();
  arma::mat t_Y_T_N = Y_N_T.t();
  arma::mat t_w_T_N = w_N_T.t();
  for(int i = 0; i < N; i++) {
    arma::vec phi_i = t_phi_3plusr_N.col(i);
    arma::vec y_i_T = t_Y_T_N.col(i);
    arma::vec w_i_T = t_w_T_N.col(i);
    double loss_i = fc_loss_QR(tau, phi_i, y_i_T, Factors_r_T, w_i_T);
    loss_PanelQR = loss_PanelQR + loss_i;
  }
  return loss_PanelQR;
}

double fc3_BIC_delta_l_tau(int l, double delta_l, double tau, arma::mat Y_N_T, arma::mat Factors_r_T, arma::mat w_N_T, arma::mat phitilde_N_3plusr, arma::mat rank_mat) {
  // determine the threshold \delta_\ell using BIC, for \ell = 1,2,3, ..., 3+r
  Environment OPT = Environment::namespace_env("PanelQGARCHwithFactor");
  Function f = OPT["fr3_BIC_delta_l_tau"];
  // Function f("fr3_BIC_delta_l_tau");
  NumericVector BIC_nv = f(l, delta_l, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat);
  arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
  double BIC_delta_l = BIC(0);
  return BIC_delta_l;
}

double fc3_select_delta_l_tau_BIC(int l, double tolerance, int num_candid_onesplit, arma::vec delta_123r_candid_lower, arma::vec delta_123r_candid_upper, double tau, arma::mat Y_N_T, arma::mat Factors_r_T, arma::mat w_N_T, arma::mat phitilde_N_3plusr, arma::mat rank_mat) {
  // select the threshold \delta_\ell using BIC, for \ell = 1,2,3, ..., 3+r
  double deltahat_l = 0.0;
  int i = 0;
  double delta_l_candid_lower = delta_123r_candid_lower(l-1);
  double delta_l_candid_upper = delta_123r_candid_upper(l-1);
  arma::vec delta_l_candids = arma::linspace(delta_l_candid_lower, delta_l_candid_upper, num_candid_onesplit);
  arma::Col<int> indexes_minBICs_from0(1);
  double minBIC_thisloop = 0.0; double minBIC_lastloop = 0.0;
  Rcout << "For \ell = " << l << ": \n";
  do {
    minBIC_lastloop = minBIC_thisloop;
    Rcout << "i = " << i << ": \n";
    Rcout << "candidates of delta_l: \n" << delta_l_candids << ". \n";
    //// minimize BIC in this loop
    arma::mat BICs_candidates(num_candid_onesplit, 1); BICs_candidates.fill(0.0);
    indexes_minBICs_from0 = {0};
    double delta_l_candid_j = 0.0;
    double BIC_before_j = 100000000.00;
    for(int j = 0; j < num_candid_onesplit; j++) {
      delta_l_candid_j = delta_l_candids(j);
      double BIC_j = fc3_BIC_delta_l_tau(l, delta_l_candid_j, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat);
      BICs_candidates(j,0) = BIC_j;
      //// search the smallest BIC and its index
      arma::Col<int> index_j = {j};
      if(BIC_j == BIC_before_j) {
        indexes_minBICs_from0 = join_cols(indexes_minBICs_from0, index_j);
      }
      if(BIC_j < BIC_before_j) {
        indexes_minBICs_from0 = index_j;
        //// update for next loop
        BIC_before_j = BIC_j;
      }
    }
    int index_first_minBIC_from0 = min(indexes_minBICs_from0); //// the index of the first smallest BICs
    int index_last_minBIC_from0 = max(indexes_minBICs_from0); //// the index of the last smallest BICs
    Rcout << "The smallest index of the smallest BICs is " << index_first_minBIC_from0 << "; ";
    Rcout << "the largest index of the smallest BICs is " << index_last_minBIC_from0 << ". \n";
    deltahat_l = delta_l_candids(index_first_minBIC_from0);
    // deltahat_l = delta_l_candids(index_last_minBIC_from0);
    Rcout << "Selected threshold deltahat_" << l << " = " << deltahat_l << ". \n";
    //// judgement for stopping this loop
    minBIC_thisloop = min(vectorise(BICs_candidates)); //// the smallest BIC 
    Rcout << "The smallest BIC is " << minBIC_thisloop << ". \n";
    //// update for next loop
    i++;
    if(index_first_minBIC_from0 == 0 && index_last_minBIC_from0 == num_candid_onesplit-1) {
      break;
    } else if(index_first_minBIC_from0 == 0) {
      delta_l_candid_lower = delta_l_candids(0);
      delta_l_candid_upper = delta_l_candids(index_last_minBIC_from0+1);
      delta_l_candids = arma::linspace(delta_l_candid_lower, delta_l_candid_upper, num_candid_onesplit+1).subvec(0, num_candid_onesplit-1);
    } else if(index_last_minBIC_from0 == num_candid_onesplit-1) {
      delta_l_candid_lower = delta_l_candids(index_first_minBIC_from0-1);
      delta_l_candid_upper = delta_l_candids(num_candid_onesplit-1);
      delta_l_candids = arma::linspace(delta_l_candid_lower, delta_l_candid_upper, num_candid_onesplit+1).subvec(1, num_candid_onesplit);
    } else {
      delta_l_candid_lower = delta_l_candids(index_first_minBIC_from0-1);
      delta_l_candid_upper = delta_l_candids(index_last_minBIC_from0+1);
      delta_l_candids = arma::linspace(delta_l_candid_lower, delta_l_candid_upper, num_candid_onesplit+2).subvec(1, num_candid_onesplit);
    }
  } while (std::abs(minBIC_thisloop - minBIC_lastloop) > tolerance);
  Rcout << "Final selected threshold deltahat_" << l << "=" << deltahat_l << ". \n";
  return deltahat_l;
}

// [[Rcpp::export]]
arma::vec fc3_1_select_delta_tau_BIC(double tolerance, int num_candid_onesplit, arma::vec delta_123r_candid_lower, arma::vec delta_123r_candid_upper, double tau, arma::mat Y_N_T, arma::mat Factors_r_T, arma::mat w_N_T, arma::mat phitilde_N_3plusr, arma::mat rank_mat) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3, ..., , \delta_{3+r})' using BIC
  int rplus3 = phitilde_N_3plusr.n_cols; 
  arma::vec deltahat_123r(rplus3); deltahat_123r.fill(0.0);
  for(int l = 1; l < (rplus3+1); l++) {
    deltahat_123r(l-1) = fc3_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123r_candid_lower, delta_123r_candid_upper, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat);
  }
  return deltahat_123r;
}

// [[Rcpp::export]]
arma::vec fc3_2_select_delta_tau_BIC(arma::Col<int> if_group_numeric, int if_mean_numeric, double tolerance, int num_candid_onesplit, arma::vec delta_123r_candid_lower, arma::vec delta_123r_candid_upper, double tau, arma::mat Y_N_T, arma::mat Factors_r_T, arma::mat w_N_T, arma::mat phitilde_N_3plusr, arma::mat rank_mat) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3, ..., , \delta_{3+r})' using BIC
  int rplus3 = phitilde_N_3plusr.n_cols; 
  arma::vec deltahat_123r(rplus3); deltahat_123r.fill(100.0);
  if(if_mean_numeric == 0) {
    for(int l = 1; l < (rplus3+1); l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123r(l-1) = fc3_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123r_candid_lower, delta_123r_candid_upper, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat);
      }
    }
  } else if(if_mean_numeric == 1) {
    int N = phitilde_N_3plusr.n_rows;
    arma::mat phitilde_N_3plusr_givenIFgroup = phitilde_N_3plusr;
    for(int l = 1; l < (rplus3+1); l++) {
      if(if_group_numeric(l-1) == 0) {
        double mean_col_l = mean(phitilde_N_3plusr.col(l-1));
        arma::mat midmat(N, 1); midmat.fill(mean_col_l); 
        phitilde_N_3plusr_givenIFgroup.col(l-1) = midmat.col(0);
      }
    }
    for(int l = 1; l < (rplus3+1); l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123r(l-1) = fc3_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123r_candid_lower, delta_123r_candid_upper, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr_givenIFgroup, rank_mat);
      }
    }
  }
  return deltahat_123r;
}

