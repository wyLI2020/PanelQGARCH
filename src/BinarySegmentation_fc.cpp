#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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

double fc_loss_QR(double tau, arma::vec phi_i, arma::vec y_i_T, arma::vec w_i_T) {
  // the loss function of the initial estimator for individual i
  Environment OPT = Environment::namespace_env("PanelQGARCH");
  Function f = OPT["Loss_QR"];
  NumericVector loss_i_nv = f(phi_i, y_i_T, w_i_T, tau);
  arma::vec loss_i(loss_i_nv.begin(), loss_i_nv.size(), false);
  return loss_i(0);
}

// [[Rcpp::export]]
double fc_loss_panelQR(double tau, arma::mat phi_N_3, arma::mat Y_N_T, arma::mat w_N_T) {
  // the loss function of the final estimator 
  int N = Y_N_T.n_rows;
  double loss_PanelQR = 0.0;
  arma::mat t_phi_3_N = phi_N_3.t();
  arma::mat t_Y_T_N = Y_N_T.t();
  arma::mat t_w_T_N = w_N_T.t();
  for(int i = 0; i < N; i++) {
    arma::vec phi_i = t_phi_3_N.col(i);
    arma::vec y_i_T = t_Y_T_N.col(i);
    arma::vec w_i_T = t_w_T_N.col(i);
    double loss_i = fc_loss_QR(tau, phi_i, y_i_T, w_i_T);
    loss_PanelQR = loss_PanelQR + loss_i;
  }
  return loss_PanelQR;
}

double fc_BIC_delta_tau(arma::vec delta_123, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // determine the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using BIC
  Environment OPT = Environment::namespace_env("PanelQGARCH");
  Function f = OPT["fr_BIC_delta_tau"];
  NumericVector BIC_nv = f(delta_123, tau, Y_N_T, w_N_T, phitilde_N_3);
  arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
  return BIC(0);
}

double fc_CV_delta_tau(arma::vec delta_123, double tau, arma::mat Y_N_T, int testLength, arma::mat w_N_T_train, arma::mat phitilde_N_3_train) {
  // determine the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' via cross-validation
  Environment OPT = Environment::namespace_env("PanelQGARCH");
  Function f = OPT["fr_CV_delta_tau"];
  NumericVector prediction_error_nv = f(delta_123, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train);
  arma::vec prediction_error(prediction_error_nv.begin(), prediction_error_nv.size(), false);
  return prediction_error(0);
}

double fc1_BIC_delta_l_tau(int l, double delta_l, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // determine the threshold \delta_\ell using BIC, for \ell = 1,2,3
  Environment OPT = Environment::namespace_env("PanelQGARCH");
  Function f1 = OPT["fr1_BIC_delta_1_tau"];
  Function f2 = OPT["fr1_BIC_delta_2_tau"];
  Function f3 = OPT["fr1_BIC_delta_3_tau"];
  double BIC_delta_l = 0.0; 
  if(l == 1) {
    NumericVector BIC_nv = f1(delta_l, tau, Y_N_T, w_N_T, phitilde_N_3);
    arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
    BIC_delta_l = BIC(0);
  } else if(l == 2) {
    NumericVector BIC_nv = f2(delta_l, tau, Y_N_T, w_N_T, phitilde_N_3);
    arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
    BIC_delta_l = BIC(0);
  } else if(l == 3) {
    NumericVector BIC_nv = f3(delta_l, tau, Y_N_T, w_N_T, phitilde_N_3);
    arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
    BIC_delta_l = BIC(0);
  }
  return BIC_delta_l;
}

double fc2_BIC_delta_l_tau(int l, double delta_l, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // determine the threshold \delta_\ell using BIC, for \ell = 1,2,3
  Environment OPT = Environment::namespace_env("PanelQGARCH");
  Function f1 = OPT["fr2_BIC_delta_1_tau"];
  Function f2 = OPT["fr2_BIC_delta_2_tau"];
  Function f3 = OPT["fr2_BIC_delta_3_tau"];
  double BIC_delta_l = 0.0; 
  if(l == 1) {
    NumericVector BIC_nv = f1(delta_l, tau, Y_N_T, w_N_T, phitilde_N_3);
    arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
    BIC_delta_l = BIC(0);
  } else if(l == 2) {
    NumericVector BIC_nv = f2(delta_l, tau, Y_N_T, w_N_T, phitilde_N_3);
    arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
    BIC_delta_l = BIC(0);
  } else if(l == 3) {
    NumericVector BIC_nv = f3(delta_l, tau, Y_N_T, w_N_T, phitilde_N_3);
    arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
    BIC_delta_l = BIC(0);
  }
  return BIC_delta_l;
}

double fc3_BIC_delta_l_tau(int l, double delta_l, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // determine the threshold \delta_\ell using BIC, for \ell = 1,2,3
  Environment OPT = Environment::namespace_env("PanelQGARCH");
  Function f1 = OPT["fr3_BIC_delta_1_tau"];
  Function f2 = OPT["fr3_BIC_delta_2_tau"];
  Function f3 = OPT["fr3_BIC_delta_3_tau"];
  double BIC_delta_l = 0.0; 
  if(l == 1) {
    NumericVector BIC_nv = f1(delta_l, tau, Y_N_T, w_N_T, phitilde_N_3);
    arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
    BIC_delta_l = BIC(0);
  } else if(l == 2) {
    NumericVector BIC_nv = f2(delta_l, tau, Y_N_T, w_N_T, phitilde_N_3);
    arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
    BIC_delta_l = BIC(0);
  } else if(l == 3) {
    NumericVector BIC_nv = f3(delta_l, tau, Y_N_T, w_N_T, phitilde_N_3);
    arma::vec BIC(BIC_nv.begin(), BIC_nv.size(), false);
    BIC_delta_l = BIC(0);
  }
  return BIC_delta_l;
}

double fc_CV_delta_l_tau(int l, double delta_l, double tau, arma::mat Y_N_T, int testLength, arma::mat w_N_T_train, arma::mat phitilde_N_3_train) {
  // determine the threshold \delta_\ell using CV, for \ell = 1,2,3
  Environment OPT = Environment::namespace_env("PanelQGARCH");
  Function f1 = OPT["fr_CV_delta_1_tau"];
  Function f2 = OPT["fr_CV_delta_2_tau"];
  Function f3 = OPT["fr_CV_delta_3_tau"];
  double CV_delta_l = 0.0; 
  if(l == 1) {
    NumericVector prediction_error_nv = f1(delta_l, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train);
    arma::vec prediction_error(prediction_error_nv.begin(), prediction_error_nv.size(), false);
    CV_delta_l = prediction_error(0);
  } else if(l == 2) {
    NumericVector prediction_error_nv = f2(delta_l, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train);
    arma::vec prediction_error(prediction_error_nv.begin(), prediction_error_nv.size(), false);
    CV_delta_l = prediction_error(0);
  } else if(l == 3) {
    NumericVector prediction_error_nv = f3(delta_l, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train);
    arma::vec prediction_error(prediction_error_nv.begin(), prediction_error_nv.size(), false);
    CV_delta_l = prediction_error(0);
  }
  return CV_delta_l;
}

// [[Rcpp::export]]
arma::vec fc_select_delta_tau_BIC(double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using BIC
  arma::vec deltahat_123(3); deltahat_123.fill(0.0);
  // int i = 0;
  // arma::vec delta_1_candids = arma::linspace(delta_123_candid_lower(0), delta_123_candid_upper(0), num_candid_onesplit+2).subvec(1, num_candid_onesplit);
  // arma::vec delta_2_candids = arma::linspace(delta_123_candid_lower(1), delta_123_candid_upper(1), num_candid_onesplit+2).subvec(1, num_candid_onesplit);
  // arma::vec delta_3_candids = arma::linspace(delta_123_candid_lower(2), delta_123_candid_upper(2), num_candid_onesplit+2).subvec(1, num_candid_onesplit);
  arma::vec delta_1_candids = arma::linspace(delta_123_candid_lower(0), delta_123_candid_upper(0), num_candid_onesplit);
  arma::vec delta_2_candids = arma::linspace(delta_123_candid_lower(1), delta_123_candid_upper(1), num_candid_onesplit);
  arma::vec delta_3_candids = arma::linspace(delta_123_candid_lower(2), delta_123_candid_upper(2), num_candid_onesplit);
  arma::Mat<int> cube_indexes_min_BICs_from0(3, 1); cube_indexes_min_BICs_from0.fill(0);
  double min_BIC_thisloop = 0.0; double min_BIC_lastloop = 0.0;
  do {
    min_BIC_lastloop = min_BIC_thisloop;
    // Rcout << "i = " << i << ": \n";
    // Rcout << "candidates of delta_1: \n" << delta_1_candids << ";";
    // Rcout << "candidates of delta_2: \n" << delta_2_candids << ";";
    // Rcout << "candidates of delta_3: \n" << delta_3_candids << ". \n";
    //// minimize BIC in this loop
    arma::vec delta_123_candid_j1j2j3(3); delta_123_candid_j1j2j3(0) = 0.0; delta_123_candid_j1j2j3(1) = 0.0; delta_123_candid_j1j2j3(2) = 0.0;
    arma::cube BICs_cube(num_candid_onesplit, num_candid_onesplit, num_candid_onesplit); BICs_cube.fill(0.0);
    // arma::Col<int> cube_index_first_min_BIC_from0(3); cube_index_first_min_BIC_from0(0) = 0; cube_index_first_min_BIC_from0(1) = 0; cube_index_first_min_BIC_from0(2) = 0; //// the index of the first smallest BICs
    // arma::Col<int> cube_index_last_min_BIC_from0(3); cube_index_last_min_BIC_from0(0) = 0; cube_index_last_min_BIC_from0(1) = 0; cube_index_last_min_BIC_from0(2) = 0; //// the index of the last smallest BICs
    arma::Col<int> zeroes_3 = {0, 0, 0}; cube_indexes_min_BICs_from0 = zeroes_3;
    double BIC_before_j1j2j3 = 100000000.00;
    for(int j1 = 0; j1 < num_candid_onesplit; j1++) {
      delta_123_candid_j1j2j3(0) = delta_1_candids(j1);
      for(int j2 = 0; j2 < num_candid_onesplit; j2++) {
        delta_123_candid_j1j2j3(1) = delta_2_candids(j2);
        for(int j3 = 0; j3 < num_candid_onesplit; j3++) {
          delta_123_candid_j1j2j3(2) = delta_3_candids(j3);
          double BIC_j1j2j3 = fc_BIC_delta_tau(delta_123_candid_j1j2j3, tau, Y_N_T, w_N_T, phitilde_N_3);
          BICs_cube(j1, j2, j3) = BIC_j1j2j3;
          //// search the smallest BIC and its index
          arma::Col<int> cube_index_j1j2j3(3); cube_index_j1j2j3(0) = j1; cube_index_j1j2j3(1) = j2; cube_index_j1j2j3(2) = j3;
          if(BIC_j1j2j3 == BIC_before_j1j2j3) {
            // cube_index_last_min_BIC_from0(0) = j1; cube_index_last_min_BIC_from0(1) = j2; cube_index_last_min_BIC_from0(2) = j3; //// the index of the last smallest BICs
            cube_indexes_min_BICs_from0 = join_rows(cube_indexes_min_BICs_from0, cube_index_j1j2j3);
          }
          if(BIC_j1j2j3 < BIC_before_j1j2j3) {
            // cube_index_first_min_BIC_from0(0) = j1; cube_index_first_min_BIC_from0(1) = j2; cube_index_first_min_BIC_from0(2) = j3; //// the index of the first smallest BICs
            cube_indexes_min_BICs_from0 = cube_index_j1j2j3;
            //// update for next loop
            BIC_before_j1j2j3 = BIC_j1j2j3;
          }
        }
      }
    }
    // if(cube_index_last_min_BIC_from0(0) == 0 && cube_index_last_min_BIC_from0(1) == 0 && cube_index_last_min_BIC_from0(2) == 0) {
    //   cube_index_last_min_BIC_from0 = cube_index_first_min_BIC_from0;
    // }
    // Rcout << "The index of the first smallest BICs in cube is \n" << cube_index_first_min_BIC_from0 << ", \n";
    // Rcout << "the index of the last smallest BICs in cube is \n" << cube_index_last_min_BIC_from0 << ". \n";
    arma::Col<int> cube_index_first_min_BIC_from0 = min(cube_indexes_min_BICs_from0, 1); //// the index of the first smallest BICs
    arma::Col<int> cube_index_last_min_BIC_from0 = max(cube_indexes_min_BICs_from0, 1); //// the index of the last smallest BICs
    // Rcout << "The indexes of the smallest BICs in cube is \n" << cube_indexes_min_BICs_from0 << ". \n";
    // Rcout << "The smallest index of the smallest BICs in cube is \n" << cube_index_first_min_BIC_from0 << "; \n";
    // Rcout << "the largest index of the smallest BICs in cube is \n" << cube_index_last_min_BIC_from0 << ". \n";
    deltahat_123(0) = delta_1_candids(cube_index_first_min_BIC_from0(0)); deltahat_123(1) = delta_2_candids(cube_index_first_min_BIC_from0(1)); deltahat_123(2) = delta_3_candids(cube_index_first_min_BIC_from0(2));
    // Rcout << "Selected thresholds: " << "deltahat_1 = " << deltahat_123(0) << "," << "deltahat_2 = " << deltahat_123(1) << "," << "deltahat_3 = " << deltahat_123(2) << ". \n";
    //// judgement for stopping this loop
    min_BIC_thisloop = min(vectorise(BICs_cube)); //// the smallest BIC 
    // Rcout << "The smallest BIC is " << min_BIC_thisloop << ". \n";
    //// update for next loop
    // i++;
    for(int l = 0; l < 3; l++) {
      arma::vec delta_l_candids(num_candid_onesplit); 
      if(l == 0) {
        delta_l_candids = delta_1_candids;
      } else if(l == 1) {
        delta_l_candids = delta_2_candids;
      } else {
        delta_l_candids = delta_3_candids;
      }
      if(cube_index_first_min_BIC_from0(l) == 0 && cube_index_last_min_BIC_from0(l) == num_candid_onesplit-1) {
        // delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_BIC_from0(l));
        // delta_123_candid_upper(l) = delta_l_candids(cube_index_first_min_BIC_from0(l));
        // delta_l_candids = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), num_candid_onesplit);
        // if(l == 0) {
        //   delta_1_candids = delta_l_candids;
        // } else if(l == 1) {
        //   delta_2_candids = delta_l_candids;
        // } else {
        //   delta_3_candids = delta_l_candids;
        // }
      } else if(cube_index_first_min_BIC_from0(l) == 0) {
        // delta_123_candid_lower(l) = delta_l_candids(cube_index_last_min_BIC_from0(l));
        delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_BIC_from0(l));
        delta_123_candid_upper(l) = delta_l_candids(cube_index_last_min_BIC_from0(l)+1);
        delta_l_candids = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), num_candid_onesplit+1).subvec(0, num_candid_onesplit-1);
        if(l == 0) {
          delta_1_candids = delta_l_candids;
        } else if(l == 1) {
          delta_2_candids = delta_l_candids;
        } else {
          delta_3_candids = delta_l_candids;
        }
      } else if(cube_index_last_min_BIC_from0(l) == num_candid_onesplit-1) {
        delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_BIC_from0(l)-1);
        // delta_123_candid_upper(l) = delta_l_candids(cube_index_first_min_BIC_from0(l));
        delta_123_candid_upper(l) = delta_l_candids(cube_index_last_min_BIC_from0(l));
        delta_l_candids = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), num_candid_onesplit+1).subvec(1, num_candid_onesplit);
        if(l == 0) {
          delta_1_candids = delta_l_candids;
        } else if(l == 1) {
          delta_2_candids = delta_l_candids;
        } else {
          delta_3_candids = delta_l_candids;
        }
      } else {
        // delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_BIC_from0(l)-1);
        // delta_123_candid_upper(l) = delta_l_candids(cube_index_first_min_BIC_from0(l));
        // arma::vec delta_l_candids_sec1 = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), floor(num_candid_onesplit/2)+1).subvec(1, floor(num_candid_onesplit/2));
        // delta_123_candid_lower(l) = delta_l_candids(cube_index_last_min_BIC_from0(l));
        // delta_123_candid_upper(l) = delta_l_candids(cube_index_last_min_BIC_from0(l)+1);
        // arma::vec delta_l_candids_sec2 = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), ceil(num_candid_onesplit/2)+1).subvec(0, ceil(num_candid_onesplit/2)-1);
        // delta_l_candids = join_cols(delta_l_candids_sec1, delta_l_candids_sec2);
        delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_BIC_from0(l)-1);
        delta_123_candid_upper(l) = delta_l_candids(cube_index_last_min_BIC_from0(l)+1);
        delta_l_candids = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), num_candid_onesplit+2).subvec(1, num_candid_onesplit);
        if(l == 0) {
          delta_1_candids = delta_l_candids;
        } else if(l == 1) {
          delta_2_candids = delta_l_candids;
        } else {
          delta_3_candids = delta_l_candids;
        }
      }
    }
  } while (std::abs(min_BIC_thisloop - min_BIC_lastloop) > tolerance);
  // Rcout << "Final selected thresholds: " << "deltahat_123: \n" << deltahat_123;
  return deltahat_123;
}

// [[Rcpp::export]]
arma::vec fc_select_delta_tau_CV(double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, int testLength, arma::mat w_N_T_train, arma::mat phitilde_N_3_train) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' via cross-validation
  arma::vec deltahat_123(3); deltahat_123.fill(0.0);
  // int i = 0;
  // arma::vec delta_1_candids = arma::linspace(delta_123_candid_lower(0), delta_123_candid_upper(0), num_candid_onesplit+2).subvec(1, num_candid_onesplit);
  // arma::vec delta_2_candids = arma::linspace(delta_123_candid_lower(1), delta_123_candid_upper(1), num_candid_onesplit+2).subvec(1, num_candid_onesplit);
  // arma::vec delta_3_candids = arma::linspace(delta_123_candid_lower(2), delta_123_candid_upper(2), num_candid_onesplit+2).subvec(1, num_candid_onesplit);
  arma::vec delta_1_candids = arma::linspace(delta_123_candid_lower(0), delta_123_candid_upper(0), num_candid_onesplit);
  arma::vec delta_2_candids = arma::linspace(delta_123_candid_lower(1), delta_123_candid_upper(1), num_candid_onesplit);
  arma::vec delta_3_candids = arma::linspace(delta_123_candid_lower(2), delta_123_candid_upper(2), num_candid_onesplit);
  arma::Mat<int> cube_indexes_min_CVs_from0(3, 1); cube_indexes_min_CVs_from0.fill(0);
  double min_CV_thisloop = 0.0; double min_CV_lastloop = 0.0;
  do {
    min_CV_lastloop = min_CV_thisloop;
    // Rcout << "i = " << i << ": \n";
    // Rcout << "candidates of delta_1: \n" << delta_1_candids << ";";
    // Rcout << "candidates of delta_2: \n" << delta_2_candids << ";";
    // Rcout << "candidates of delta_3: \n" << delta_3_candids << ". \n";
    //// minimize CV in this loop
    arma::vec delta_123_candid_j1j2j3(3); delta_123_candid_j1j2j3(0) = 0.0; delta_123_candid_j1j2j3(1) = 0.0; delta_123_candid_j1j2j3(2) = 0.0;
    arma::cube CVs_cube(num_candid_onesplit, num_candid_onesplit, num_candid_onesplit); CVs_cube.fill(0.0);
    // arma::Col<int> cube_index_first_min_CV_from0(3); cube_index_first_min_CV_from0(0) = 0; cube_index_first_min_CV_from0(1) = 0; cube_index_first_min_CV_from0(2) = 0; //// the index of the first smallest CVs
    // arma::Col<int> cube_index_last_min_CV_from0(3); cube_index_last_min_CV_from0(0) = 0; cube_index_last_min_CV_from0(1) = 0; cube_index_last_min_CV_from0(2) = 0; //// the index of the last smallest CVs
    arma::Col<int> zeroes_3 = {0, 0, 0}; cube_indexes_min_CVs_from0 = zeroes_3;
    double CV_before_j1j2j3 = 100000000.00;
    for(int j1 = 0; j1 < num_candid_onesplit; j1++) {
      delta_123_candid_j1j2j3(0) = delta_1_candids(j1);
      for(int j2 = 0; j2 < num_candid_onesplit; j2++) {
        delta_123_candid_j1j2j3(1) = delta_2_candids(j2);
        for(int j3 = 0; j3 < num_candid_onesplit; j3++) {
          delta_123_candid_j1j2j3(2) = delta_3_candids(j3);
          double CV_j1j2j3 = fc_CV_delta_tau(delta_123_candid_j1j2j3, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train);
          CVs_cube(j1, j2, j3) = CV_j1j2j3;
          //// search the smallest CV and its index
          arma::Col<int> cube_index_j1j2j3(3); cube_index_j1j2j3(0) = j1; cube_index_j1j2j3(1) = j2; cube_index_j1j2j3(2) = j3;
          if(CV_j1j2j3 == CV_before_j1j2j3) {
            // cube_index_last_min_CV_from0(0) = j1; cube_index_last_min_CV_from0(1) = j2; cube_index_last_min_CV_from0(2) = j3; //// the index of the last smallest CVs
            cube_indexes_min_CVs_from0 = join_rows(cube_indexes_min_CVs_from0, cube_index_j1j2j3);
          }
          if(CV_j1j2j3 < CV_before_j1j2j3) {
            // cube_index_first_min_CV_from0(0) = j1; cube_index_first_min_CV_from0(1) = j2; cube_index_first_min_CV_from0(2) = j3; //// the index of the first smallest CVs
            cube_indexes_min_CVs_from0 = cube_index_j1j2j3;
            //// update for next loop
            CV_before_j1j2j3 = CV_j1j2j3;
          }
        }
      }
    }
    // if(cube_index_last_min_CV_from0(0) == 0 && cube_index_last_min_CV_from0(1) == 0 && cube_index_last_min_CV_from0(2) == 0) {
    //   cube_index_last_min_CV_from0 = cube_index_first_min_CV_from0;
    // }
    // Rcout << "The index of the first smallest CVs in cube is \n" << cube_index_first_min_CV_from0 << ", \n";
    // Rcout << "the index of the last smallest CVs in cube is \n" << cube_index_last_min_CV_from0 << ". \n";
    arma::Col<int> cube_index_first_min_CV_from0 = min(cube_indexes_min_CVs_from0, 1); //// the index of the first smallest CVs
    arma::Col<int> cube_index_last_min_CV_from0 = max(cube_indexes_min_CVs_from0, 1); //// the index of the last smallest CVs
    // Rcout << "The indexes of the smallest CVs in cube is \n" << cube_indexes_min_CVs_from0 << ". \n";
    // Rcout << "The smallest index of the smallest CVs in cube is \n" << cube_index_first_min_CV_from0 << "; \n";
    // Rcout << "the largest index of the smallest CVs in cube is \n" << cube_index_last_min_CV_from0 << ". \n";
    deltahat_123(0) = delta_1_candids(cube_index_first_min_CV_from0(0)); deltahat_123(1) = delta_2_candids(cube_index_first_min_CV_from0(1)); deltahat_123(2) = delta_3_candids(cube_index_first_min_CV_from0(2));
    // Rcout << "Selected thresholds: " << "deltahat_1 = " << deltahat_123(0) << "," << "deltahat_2 = " << deltahat_123(1) << "," << "deltahat_3 = " << deltahat_123(2) << ". \n";
    //// judgement for stopping this loop
    min_CV_thisloop = min(vectorise(CVs_cube)); //// the smallest CV 
    // Rcout << "The smallest CV is " << min_CV_thisloop << ". \n";
    //// update for next loop
    // i++;
    for(int l = 0; l < 3; l++) {
      arma::vec delta_l_candids(num_candid_onesplit); 
      if(l == 0) {
        delta_l_candids = delta_1_candids;
      } else if(l == 1) {
        delta_l_candids = delta_2_candids;
      } else {
        delta_l_candids = delta_3_candids;
      }
      if(cube_index_first_min_CV_from0(l) == 0 && cube_index_last_min_CV_from0(l) == num_candid_onesplit-1) {
        // delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_CV_from0(l));
        // delta_123_candid_upper(l) = delta_l_candids(cube_index_first_min_CV_from0(l));
        // delta_l_candids = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), num_candid_onesplit);
        // if(l == 0) {
        //   delta_1_candids = delta_l_candids;
        // } else if(l == 1) {
        //   delta_2_candids = delta_l_candids;
        // } else {
        //   delta_3_candids = delta_l_candids;
        // }
      } else if(cube_index_first_min_CV_from0(l) == 0) {
        // delta_123_candid_lower(l) = delta_l_candids(cube_index_last_min_CV_from0(l));
        delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_CV_from0(l));
        delta_123_candid_upper(l) = delta_l_candids(cube_index_last_min_CV_from0(l)+1);
        delta_l_candids = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), num_candid_onesplit+1).subvec(0, num_candid_onesplit-1);
        if(l == 0) {
          delta_1_candids = delta_l_candids;
        } else if(l == 1) {
          delta_2_candids = delta_l_candids;
        } else {
          delta_3_candids = delta_l_candids;
        }
      } else if(cube_index_last_min_CV_from0(l) == num_candid_onesplit-1) {
        delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_CV_from0(l)-1);
        // delta_123_candid_upper(l) = delta_l_candids(cube_index_first_min_CV_from0(l));
        delta_123_candid_upper(l) = delta_l_candids(cube_index_last_min_CV_from0(l));
        delta_l_candids = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), num_candid_onesplit+1).subvec(1, num_candid_onesplit);
        if(l == 0) {
          delta_1_candids = delta_l_candids;
        } else if(l == 1) {
          delta_2_candids = delta_l_candids;
        } else {
          delta_3_candids = delta_l_candids;
        }
      } else {
        // delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_CV_from0(l)-1);
        // delta_123_candid_upper(l) = delta_l_candids(cube_index_first_min_CV_from0(l));
        // arma::vec delta_l_candids_sec1 = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), floor(num_candid_onesplit/2)+1).subvec(1, floor(num_candid_onesplit/2));
        // delta_123_candid_lower(l) = delta_l_candids(cube_index_last_min_CV_from0(l));
        // delta_123_candid_upper(l) = delta_l_candids(cube_index_last_min_CV_from0(l)+1);
        // arma::vec delta_l_candids_sec2 = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), ceil(num_candid_onesplit/2)+1).subvec(0, ceil(num_candid_onesplit/2)-1);
        // delta_l_candids = join_cols(delta_l_candids_sec1, delta_l_candids_sec2);
        delta_123_candid_lower(l) = delta_l_candids(cube_index_first_min_CV_from0(l)-1);
        delta_123_candid_upper(l) = delta_l_candids(cube_index_last_min_CV_from0(l)+1);
        delta_l_candids = arma::linspace(delta_123_candid_lower(l), delta_123_candid_upper(l), num_candid_onesplit+2).subvec(1, num_candid_onesplit);
        if(l == 0) {
          delta_1_candids = delta_l_candids;
        } else if(l == 1) {
          delta_2_candids = delta_l_candids;
        } else {
          delta_3_candids = delta_l_candids;
        }
      }
    }
  } while (std::abs(min_CV_thisloop - min_CV_lastloop) > tolerance);
  // Rcout << "Final selected thresholds: " << "deltahat_123: \n" << deltahat_123;
  return deltahat_123;
}

double fc1_select_delta_l_tau_BIC(int l, double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \delta_\ell using BIC, for \ell = 1,2,3
  double deltahat_l = 0.0;
  int i = 0;
  double delta_l_candid_lower = delta_123_candid_lower(l-1);
  double delta_l_candid_upper = delta_123_candid_upper(l-1);
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
      double BIC_j = fc1_BIC_delta_l_tau(l, delta_l_candid_j, tau, Y_N_T, w_N_T, phitilde_N_3);
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
arma::vec fc1_1_select_delta_tau_BIC(double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using BIC
  arma::vec deltahat_123(3); deltahat_123.fill(0.0);
  for(int l = 1; l < 4; l++) {
    deltahat_123(l-1) = fc1_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3);
  }
  return deltahat_123;
}

// [[Rcpp::export]]
arma::vec fc1_2_select_delta_tau_BIC(arma::Col<int> if_group_numeric, int if_mean_numeric, double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using BIC
  arma::vec deltahat_123(3); deltahat_123.fill(100.0);
  if(if_mean_numeric == 0) {
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123(l-1) = fc1_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3);
      }
    }
  } else if(if_mean_numeric == 1) {
    int N = phitilde_N_3.n_rows;
    arma::mat phitilde_N_3_givenIFgroup = phitilde_N_3;
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 0) {
        double mean_col_l = mean(phitilde_N_3.col(l-1));
        arma::mat midmat(N, 1); midmat.fill(mean_col_l); 
        phitilde_N_3_givenIFgroup.col(l-1) = midmat.col(0);
      }
    }
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123(l-1) = fc1_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3_givenIFgroup);
      }
    }
  }
  return deltahat_123;
}

double fc2_select_delta_l_tau_BIC(int l, double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \delta_\ell using BIC, for \ell = 1,2,3
  double deltahat_l = 0.0;
  int i = 0;
  double delta_l_candid_lower = delta_123_candid_lower(l-1);
  double delta_l_candid_upper = delta_123_candid_upper(l-1);
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
      double BIC_j = fc2_BIC_delta_l_tau(l, delta_l_candid_j, tau, Y_N_T, w_N_T, phitilde_N_3);
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
arma::vec fc2_1_select_delta_tau_BIC(double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using BIC
  arma::vec deltahat_123(3); deltahat_123.fill(0.0);
  for(int l = 1; l < 4; l++) {
    deltahat_123(l-1) = fc2_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3);
  }
  return deltahat_123;
}

// [[Rcpp::export]]
arma::vec fc2_2_select_delta_tau_BIC(arma::Col<int> if_group_numeric, int if_mean_numeric, double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using BIC
  arma::vec deltahat_123(3); deltahat_123.fill(100.0);
  if(if_mean_numeric == 0) {
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123(l-1) = fc2_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3);
      }
    }
  } else if(if_mean_numeric == 1) {
    int N = phitilde_N_3.n_rows;
    arma::mat phitilde_N_3_givenIFgroup = phitilde_N_3;
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 0) {
        double mean_col_l = mean(phitilde_N_3.col(l-1));
        arma::mat midmat(N, 1); midmat.fill(mean_col_l); 
        phitilde_N_3_givenIFgroup.col(l-1) = midmat.col(0);
      }
    }
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123(l-1) = fc2_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3_givenIFgroup);
      }
    }
  }
  return deltahat_123;
}

double fc3_select_delta_l_tau_BIC(int l, double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \delta_\ell using BIC, for \ell = 1,2,3
  double deltahat_l = 0.0;
  int i = 0;
  double delta_l_candid_lower = delta_123_candid_lower(l-1);
  double delta_l_candid_upper = delta_123_candid_upper(l-1);
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
      double BIC_j = fc3_BIC_delta_l_tau(l, delta_l_candid_j, tau, Y_N_T, w_N_T, phitilde_N_3);
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
arma::vec fc3_1_select_delta_tau_BIC(double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using BIC
  arma::vec deltahat_123(3); deltahat_123.fill(0.0);
  for(int l = 1; l < 4; l++) {
    deltahat_123(l-1) = fc3_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3);
  }
  return deltahat_123;
}

// [[Rcpp::export]]
arma::vec fc3_2_select_delta_tau_BIC(arma::Col<int> if_group_numeric, int if_mean_numeric, double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, arma::mat w_N_T, arma::mat phitilde_N_3) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using BIC
  arma::vec deltahat_123(3); deltahat_123.fill(100.0);
  if(if_mean_numeric == 0) {
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123(l-1) = fc3_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3);
      }
    }
  } else if(if_mean_numeric == 1) {
    int N = phitilde_N_3.n_rows;
    arma::mat phitilde_N_3_givenIFgroup = phitilde_N_3;
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 0) {
        double mean_col_l = mean(phitilde_N_3.col(l-1));
        arma::mat midmat(N, 1); midmat.fill(mean_col_l); 
        phitilde_N_3_givenIFgroup.col(l-1) = midmat.col(0);
      }
    }
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123(l-1) = fc3_select_delta_l_tau_BIC(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3_givenIFgroup);
      }
    }
  }
  return deltahat_123;
}

double fc_select_delta_l_tau_CV(int l, double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, int testLength, arma::mat w_N_T_train, arma::mat phitilde_N_3_train) {
  // select the threshold \delta_\ell using CV, for \ell = 1,2,3
  double deltahat_l = 0.0;
  int i = 0;
  double delta_l_candid_lower = delta_123_candid_lower(l-1);
  double delta_l_candid_upper = delta_123_candid_upper(l-1);
  arma::vec delta_l_candids = arma::linspace(delta_l_candid_lower, delta_l_candid_upper, num_candid_onesplit);
  arma::Col<int> indexes_minCVs_from0(1);
  double minCV_thisloop = 0.0; double minCV_lastloop = 0.0;
  Rcout << "For \ell = " << l << ": \n";
  do {
    minCV_lastloop = minCV_thisloop;
    Rcout << "i = " << i << ": \n";
    Rcout << "candidates of delta_l: \n" << delta_l_candids << ". \n";
    //// minimize CV in this loop
    arma::mat CVs_candidates(num_candid_onesplit, 1); CVs_candidates.fill(0.0);
    indexes_minCVs_from0 = {0};
    double delta_l_candid_j = 0.0;
    double CV_before_j = 100000000.00;
    for(int j = 0; j < num_candid_onesplit; j++) {
      delta_l_candid_j = delta_l_candids(j);
      double CV_j = fc_CV_delta_l_tau(l, delta_l_candid_j, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train);
      CVs_candidates(j,0) = CV_j;
      //// search the smallest CV and its index
      arma::Col<int> index_j = {j};
      if(CV_j == CV_before_j) {
        indexes_minCVs_from0 = join_cols(indexes_minCVs_from0, index_j);
      }
      if(CV_j < CV_before_j) {
        indexes_minCVs_from0 = index_j;
        //// update for next loop
        CV_before_j = CV_j;
      }
    }
    int index_first_minCV_from0 = min(indexes_minCVs_from0); //// the index of the first smallest CVs
    int index_last_minCV_from0 = max(indexes_minCVs_from0); //// the index of the last smallest CVs
    Rcout << "The smallest index of the smallest CVs is " << index_first_minCV_from0 << "; ";
    Rcout << "the largest index of the smallest CVs is " << index_last_minCV_from0 << ". \n";
    deltahat_l = delta_l_candids(index_first_minCV_from0);
    // deltahat_l = delta_l_candids(index_last_minCV_from0);
    Rcout << "Selected threshold deltahat_" << l << " = " << deltahat_l << ". \n";
    //// judgement for stopping this loop
    minCV_thisloop = min(vectorise(CVs_candidates)); //// the smallest CV 
    Rcout << "The smallest CV is " << minCV_thisloop << ". \n";
    //// update for next loop
    i++;
    if(index_first_minCV_from0 == 0 && index_last_minCV_from0 == num_candid_onesplit-1) {
      break;
    } else if(index_first_minCV_from0 == 0) {
      delta_l_candid_lower = delta_l_candids(0);
      delta_l_candid_upper = delta_l_candids(index_last_minCV_from0+1);
      delta_l_candids = arma::linspace(delta_l_candid_lower, delta_l_candid_upper, num_candid_onesplit+1).subvec(0, num_candid_onesplit-1);
    } else if(index_last_minCV_from0 == num_candid_onesplit-1) {
      delta_l_candid_lower = delta_l_candids(index_first_minCV_from0-1);
      delta_l_candid_upper = delta_l_candids(num_candid_onesplit-1);
      delta_l_candids = arma::linspace(delta_l_candid_lower, delta_l_candid_upper, num_candid_onesplit+1).subvec(1, num_candid_onesplit);
    } else {
      delta_l_candid_lower = delta_l_candids(index_first_minCV_from0-1);
      delta_l_candid_upper = delta_l_candids(index_last_minCV_from0+1);
      delta_l_candids = arma::linspace(delta_l_candid_lower, delta_l_candid_upper, num_candid_onesplit+2).subvec(1, num_candid_onesplit);
    }
  } while (std::abs(minCV_thisloop - minCV_lastloop) > tolerance);
  Rcout << "Final selected threshold deltahat_" << l << "=" << deltahat_l << ". \n";
  return deltahat_l;
}

// [[Rcpp::export]]
arma::vec fc_1_select_delta_tau_CV(double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, int testLength, arma::mat w_N_T_train, arma::mat phitilde_N_3_train) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using CV
  arma::vec deltahat_123(3); deltahat_123.fill(0.0);
  for(int l = 1; l < 4; l++) {
    deltahat_123(l-1) = fc_select_delta_l_tau_CV(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train);
  }
  return deltahat_123;
}

// [[Rcpp::export]]
arma::vec fc_2_select_delta_tau_CV(arma::Col<int> if_group_numeric, int if_mean_numeric, double tolerance, int num_candid_onesplit, arma::vec delta_123_candid_lower, arma::vec delta_123_candid_upper, double tau, arma::mat Y_N_T, int testLength, arma::mat w_N_T_train, arma::mat phitilde_N_3_train) {
  // select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using CV
  arma::vec deltahat_123(3); deltahat_123.fill(100.0);
  if(if_mean_numeric == 0) {
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123(l-1) = fc_select_delta_l_tau_CV(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train);
      }
    }
  } else if(if_mean_numeric == 1) {
    int N = phitilde_N_3_train.n_rows;
    arma::mat phitilde_N_3_train_givenIFgroup = phitilde_N_3_train;
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 0) {
        double mean_col_l = mean(phitilde_N_3_train.col(l-1));
        arma::mat midmat(N, 1); midmat.fill(mean_col_l); 
        phitilde_N_3_train_givenIFgroup.col(l-1) = midmat.col(0);
      }
    }
    for(int l = 1; l < 4; l++) {
      if(if_group_numeric(l-1) == 1) {
        deltahat_123(l-1) = fc_select_delta_l_tau_CV(l, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train_givenIFgroup);
      }
    }
  }
  return deltahat_123;
}

