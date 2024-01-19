# Estimation Procedure

## Initial estimates
# library(rugarch)
# library(dfoptim)
# source("InitialEstimation_fr.R")
# Rcpp::sourceCpp("InitialEstimation_fc.cpp")

fr_weights_NT <- function(Y_N_T) {
  # matrix of weights, i.e. [w_{it}] for i=1,...,N and t=1,...,T
  c_N <- apply(Y_N_T, MARGIN=1, FUN=quantile, probs=c(0.95))
  w_N_T <- fc_weighs_NT(Y_N_T, c_N)
  return(w_N_T)
}


## Change points detection by binary segmentation algorithm
# source("BinarySegmentation_fr.R")
# Rcpp::sourceCpp("BinarySegmentation_fc.cpp")

fr_select_delta_tau <- function(if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3, testLength, w_N_T_train, phitilde_N_3_train) {
  # select the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)'
  if(is.na(if_group[1])) {
    if(criterion == "BIC1") {
      deltahat_123 <- fc1_1_select_delta_tau_BIC(tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3)
      return(deltahat_123)
    } else if(criterion == "BIC2") {
      deltahat_123 <- fc2_1_select_delta_tau_BIC(tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3)
      return(deltahat_123)
    } else if(criterion == "BIC3") {
      deltahat_123 <- fc3_1_select_delta_tau_BIC(tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3)
      return(deltahat_123)
    } else if(criterion == "CV") {
      deltahat_123 <- fc_1_select_delta_tau_CV(tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train)
      return(deltahat_123)
    }
  } else {
    if_group_numeric <- as.numeric(if_group)
    if_mean_numeric <- as.numeric(if_mean)
    if(criterion == "BIC1") {
      deltahat_123 <- fc1_2_select_delta_tau_BIC(if_group_numeric, if_mean_numeric, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3)
      return(deltahat_123)
    } else if(criterion == "BIC2") {
      deltahat_123 <- fc2_2_select_delta_tau_BIC(if_group_numeric, if_mean_numeric, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3)
      return(deltahat_123)
    } else if(criterion == "BIC3") {
      deltahat_123 <- fc3_2_select_delta_tau_BIC(if_group_numeric, if_mean_numeric, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3)
      return(deltahat_123)
    } else if(criterion == "CV") {
      deltahat_123 <- fc_2_select_delta_tau_CV(if_group_numeric, if_mean_numeric, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train)
      return(deltahat_123)
    }
  }
}

fr_partition_phi_N <- function(g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) {
  # obtain the partition of (\bm\phi_1, ..., \bm\phi_N)' from the change points and the ranks
  partition_omega_N_ascend <- rep(1:g_1, diff(c(0, changepoints_1)))
  partition_alpha_N_ascend <- rep(1:g_2, diff(c(0, changepoints_2)))
  partition_beta_N_ascend <- rep(1:g_3, diff(c(0, changepoints_3)))
  partition_omega_N <- partition_omega_N_ascend[rank_1]
  partition_alpha_N <- partition_alpha_N_ascend[rank_2]
  partition_beta_N <- partition_beta_N_ascend[rank_3]
  partition_phi_N_3 <- cbind(partition_omega_N, partition_alpha_N, partition_beta_N)
  return(partition_phi_N_3)
}

## Final estimates
fr_onerep_tau <- function(tau, Y_N_T, if_group=c(TRUE, TRUE, TRUE), if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123_candid_lower=c(0,0,0), delta_123_candid_upper=c(2,2,2), testLength=NA) {
  # {\widetilde\bm\phi}, {detected partition of \bm\phi} and {\widehat\bm\phi} obtained by the three-stage estimation procedure
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## Stage 1 - Initial estimates
  w_N_T <- fr_weights_NT(Y_N_T)
  phitilde_N_3 <- fc_phitilde_N_tau(tau, Y_N_T, w_N_T)
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## Stage 2 - Change points detection by binary segmentation algorithm
  ### threshold selection
  w_N_T_train <- NA; phitilde_N_3_train <- NA
  if(criterion == "CV") { # initial estimates (training set)
    trainLength <- T - testLength; Y_N_T_train <- Y_N_T[,1:trainLength]; Y_N_T_test <- Y_N_T[,trainLength+c(1:testLength)]
    w_N_T_train <- fr_weights_NT(Y_N_T_train)
    phitilde_N_3_train <- fc_phitilde_N_tau(tau, Y_N_T_train, w_N_T_train)
  }
  deltahat_123 <- fr_select_delta_tau(if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3, testLength, w_N_T_train, phitilde_N_3_train)
  ### change points detection
  ghat_1 <- 1; changepoints_1 <- c(N); rank_1 <- rank(omegatilde_N)
  ghat_2 <- 1; changepoints_2 <- c(N); rank_2 <- rank(alphatilde_N)
  ghat_3 <- 1; changepoints_3 <- c(N); rank_3 <- rank(betatilde_N)
  if(if_group[1] == TRUE) {
    omegatilde_N_ascend <- sort(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
    delta_omega <- deltahat_123[1]; changepoints_omegatilde <- fc_changepoints(delta_omega, 1, N, omegatilde_N_ascend); num_changepoints_omegatilde <- length(changepoints_omegatilde) # (\widehat{h}_{1,1}, ..., \widehat{h}_{1,\widehat{g}_1})' and \widehat{g}_1 
    ghat_1 <- num_changepoints_omegatilde; changepoints_1 <- changepoints_omegatilde
  }
  if(if_group[2] == TRUE) {
    alphatilde_N_ascend <- sort(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
    delta_alpha <- deltahat_123[2]; changepoints_alphatilde <- fc_changepoints(delta_alpha, 1, N, alphatilde_N_ascend); num_changepoints_alphatilde <- length(changepoints_alphatilde) # (\widehat{h}_{2,1}, ..., \widehat{h}_{2,\widehat{g}_1})' and \widehat{g}_2 
    ghat_2 <- num_changepoints_alphatilde; changepoints_2 <- changepoints_alphatilde
  }
  if(if_group[3] == TRUE) {
    betatilde_N_ascend <- sort(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
    delta_beta <- deltahat_123[3]; changepoints_betatilde <- fc_changepoints(delta_beta, 1, N, betatilde_N_ascend); num_changepoints_betatilde <- length(changepoints_betatilde) # (\widehat{h}_{3,1}, ..., \widehat{h}_{3,\widehat{g}_1})' and \widehat{g}_3 
    ghat_3 <- num_changepoints_betatilde; changepoints_3 <- changepoints_betatilde
  }
  ### partition detection
  partition_phihat_N_3 <- fr_partition_phi_N(ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  ## Stage 3 - Final estimates
  thetahat <- fr_thetahat_tau(tau, Y_N_T, w_N_T, phitilde_N_3, ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  phihat_N_3 <- fr_phi_N_homo(thetahat, ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  ## Summary
  onerep_results <- c(as.vector(phitilde_N_3), as.vector(partition_phihat_N_3), as.vector(phihat_N_3))
  names(onerep_results) <- c(paste0("omegatilde_", 1:N), paste0("alphatilde_", 1:N), paste0("betatilde_", 1:N), 
                             paste0("label_omega_", 1:N), paste0("label_alpha_", 1:N), paste0("label_beta_", 1:N), 
                             paste0("omegahat_", 1:N), paste0("alphahat_", 1:N), paste0("betahat_", 1:N))
  return(onerep_results)
}

fr_oracle_onerep_tau <- function(tau, Y_N_T, phiini_N_3, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3) {
  ## oracle estimation
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  w_N_T <- fr_weights_NT(Y_N_T)
  omegaini_N <- phiini_N_3[,1]; alphaini_N <- phiini_N_3[,2]; betaini_N <- phiini_N_3[,3]
  rank_1 <- rank(omegaini_N); rank_2 <- rank(alphaini_N); rank_3 <- rank(betaini_N)
  thetahat <- fr_thetahat_tau(tau, Y_N_T, w_N_T, phiini_N_3, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  phihat_N_3 <- fr_phi_N_homo(thetahat, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  ## Summary
  onerep_results <- as.vector(phihat_N_3)
  names(onerep_results) <- c(paste0("oracle_omegahat_", 1:N), paste0("oracle_alphahat_", 1:N), paste0("oracle_betahat_", 1:N))
  return(onerep_results)
}

fr_onerep_givenGroups_tau <- function(tau, Y_N_T, GroupLables_N_3) {
  # {\widetilde\bm\phi}, {detected partition of \bm\phi} and {\widehat\bm\phi} obtained by the three-stage estimation procedure
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## Initial estimates
  w_N_T <- fr_weights_NT(Y_N_T)
  phitilde_N_3 <- fc_phitilde_N_tau(tau, Y_N_T, w_N_T)
  lables_omega <- GroupLables_N_3[,1]
  lables_alpha <- GroupLables_N_3[,2]
  lables_beta <- GroupLables_N_3[,3]
  ## Final estimates
  g_1 <- max(lables_omega); g_2 <- max(lables_alpha); g_3 <- max(lables_beta)
  lables_omega_ascend <- sort(lables_omega); lables_alpha_ascend <- sort(lables_alpha); lables_beta_ascend <- sort(lables_beta)
  changepoints_1 <- fc_changepoints(0.0001, 1, N, lables_omega_ascend); changepoints_2 <- fc_changepoints(0.0001, 1, N, lables_alpha_ascend); changepoints_3 <- fc_changepoints(0.0001, 1, N, lables_beta_ascend); 
  rank_1 <- rank(lables_omega); rank_2 <- rank(lables_alpha); rank_3 <- rank(lables_beta)
  thetahat <- fr_thetahat_tau(tau, Y_N_T, w_N_T, phitilde_N_3, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  phihat_N_3 <- fr_phi_N_homo(thetahat, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  ## Summary
  onerep_results <- c(as.vector(phitilde_N_3), as.vector(GroupLables_N_3), as.vector(phihat_N_3))
  names(onerep_results) <- c(paste0("omegatilde_", 1:N), paste0("alphatilde_", 1:N), paste0("betatilde_", 1:N), 
                             paste0("label_omega_", 1:N), paste0("label_alpha_", 1:N), paste0("label_beta_", 1:N), 
                             paste0("omegahat_", 1:N), paste0("alphahat_", 1:N), paste0("betahat_", 1:N))
  return(onerep_results)
}


