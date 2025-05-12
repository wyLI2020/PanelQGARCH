# Estimation Procedure

## Initial estimates
# library(rugarch)
# library(dfoptim)
# source("InitialEstimation_factor_fr.R")
# Rcpp::sourceCpp("InitialEstimation_factor_fc.cpp")

fr_weights_NT <- function(Y_N_T) {
  # matrix of weights, i.e. [w_{it}] for i=1,...,N and t=1,...,T
  c_N <- apply(Y_N_T, MARGIN=1, FUN=quantile, probs=c(0.95))
  w_N_T <- fc_weighs_NT(Y_N_T, c_N)
  return(w_N_T)
}


## Change points detection by binary segmentation algorithm
# source("BinarySegmentation_factor_fr_ini1.R")
# Rcpp::sourceCpp("BinarySegmentation_factor_fc.cpp")

fr_select_delta_tau <- function(if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123r_candid_lower, delta_123r_candid_upper, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat) {
  # select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)'
  if(is.na(if_group[1])) {
    if(criterion == "BIC3") {
      deltahat_123r <- fc3_1_select_delta_tau_BIC(tolerance, num_candid_onesplit, delta_123r_candid_lower, delta_123r_candid_upper, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat)
      return(deltahat_123r)
    }
  } else {
    if_group_numeric <- as.numeric(if_group)
    if_mean_numeric <- as.numeric(if_mean)
    if(criterion == "BIC3") {
      deltahat_123r <- fc3_2_select_delta_tau_BIC(if_group_numeric, if_mean_numeric, tolerance, num_candid_onesplit, delta_123r_candid_lower, delta_123r_candid_upper, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat)
      return(deltahat_123r)
    }
  }
}

fr_partition_phi_l_N <- function(l, g_vec, changepoints_vec, rank_mat) {
  # obtain the partition of (\bm\phi_{l1}, ..., \bm\phi_{lN})' from the change points and the ranks, when l = 1,2,3, ..., 3+r
  g_l <- g_vec[l]
  if(l == 1) {
    changepoints_l <- changepoints_vec[1:g_l]
  } else {
    changepoints_l <- changepoints_vec[sum(g_vec[1:(l-1)])+(1:g_l)]
  }
  rank_l <- rank_mat[,l]
  partition_phi_l_N_ascend <- rep(1:g_l, diff(c(0, changepoints_l)))
  partition_phi_l_N <- partition_phi_l_N_ascend[rank_l]
  return(partition_phi_l_N)
}

# Rcpp::sourceCpp("EstimationProcedure_factor_fc.cpp")

## Final estimates
fr_onerep_tau <- function(tau, Y_N_T, Factors_r_T, if_group, if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123r_candid_lower, delta_123r_candid_upper) {
  # {\widetilde\bm\phi}, {detected partition of \bm\phi} and {\widehat\bm\phi} obtained by the three-stage estimation procedure
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  r <- nrow(Factors_r_T)
  ## Stage 1 - Initial estimates
  w_N_T <- fr_weights_NT(Y_N_T)
  phitilde_N_3plusr <- fc_phitilde_N_tau(tau, Y_N_T, Factors_r_T, w_N_T)
  ## Stage 2 - Change points detection by binary segmentation algorithm
  ### threshold selection
  rank_mat <- fc_rank_mat(phitilde_N_3plusr)
  deltahat_123r <- fr_select_delta_tau(if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123r_candid_lower, delta_123r_candid_upper, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat)
  ### change points detection
  ghat_changepoints_vec <- fc_g_changepoints_vec(phitilde_N_3plusr, deltahat_123r)
  ghat_vec <- ghat_changepoints_vec[1:(3+r)]
  changepoints_vec <- ghat_changepoints_vec[-c(1:(3+r))]
  ### partition detection
  partition_phihat_N_3plusr <- fc_partition_phi_N(ghat_vec, changepoints_vec, rank_mat)
  ## Stage 3 - Final estimates
  thetahat <- fr_thetahat_tau(tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, ghat_vec, changepoints_vec, rank_mat)
  phihat_N_3plusr <- fc_phi_N_homo(thetahat, ghat_vec, changepoints_vec, rank_mat)
  ## Summary
  onerep_results <- c(as.vector(phitilde_N_3plusr), as.vector(partition_phihat_N_3plusr), as.vector(phihat_N_3plusr))
  names(onerep_results) <- c(paste0("omegatilde_", 1:N), paste0("alphatilde_", 1:N), paste0("betatilde_", 1:N), paste0("lambda", rep(1:r, each=N), "tilde_", rep(1:N, times = r)), 
                             paste0("label_omega_", 1:N), paste0("label_alpha_", 1:N), paste0("label_beta_", 1:N), paste0("label_lambda", rep(1:r, each=N), "_", rep(1:N, times = r)), 
                             paste0("omegahat_", 1:N), paste0("alphahat_", 1:N), paste0("betahat_", 1:N), paste0("lambda", rep(1:r, each=N), "hat_", rep(1:N, times = r)))
  return(onerep_results)
}

