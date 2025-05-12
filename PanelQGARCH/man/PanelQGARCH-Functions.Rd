\name{PanelQGARCH-Functions}
\alias{check_function}
\alias{psi_function}
\alias{h_compute}
\alias{q_y}
\alias{g_tau}
\alias{arch_fft}
\alias{Init_par}
\alias{Loss_QR}
\alias{Loss_QR_gr}
\alias{fit1_optim}
\alias{fr_phitilde_i_tau}
\alias{fit1_optim_grid}
\alias{Sigma_QR}
\alias{ASD_QR}
\alias{fc_weighs_NT}
\alias{fc_phitilde_N_tau}
\alias{fc_q_y_t0plus1}
\alias{fr_weights_NT}
\alias{fr_phi_N_homo}
\alias{fr_loss_panelQR}
\alias{fr2_loss_panelQR}
\alias{fr_theta_initial_IN_optim}
\alias{fr_thetahat_tau}
\alias{fr_thetahat_omega_tau}
\alias{fr_thetahat_alpha_tau}
\alias{fr_thetahat_beta_tau}
\alias{fr_BIC_delta_tau}
\alias{fr_CV_delta_tau}
\alias{fr1_BIC_delta_1_tau}
\alias{fr1_BIC_delta_2_tau}
\alias{fr1_BIC_delta_3_tau}
\alias{fr2_BIC_delta_1_tau}
\alias{fr2_BIC_delta_2_tau}
\alias{fr2_BIC_delta_3_tau}
\alias{fr3_BIC_delta_1_tau}
\alias{fr3_BIC_delta_2_tau}
\alias{fr3_BIC_delta_3_tau}
\alias{fr_CV_delta_1_tau}
\alias{fr_CV_delta_2_tau}
\alias{fr_CV_delta_3_tau}
\alias{fc_calculate_segment_means}
\alias{fc_Delta_ijh}
\alias{fc_changepoints}
\alias{fc_loss_panelQR}
\alias{fc_select_delta_tau_BIC}
\alias{fc_select_delta_tau_CV}
\alias{fc1_1_select_delta_tau_BIC}
\alias{fc1_2_select_delta_tau_BIC}
\alias{fc2_1_select_delta_tau_BIC}
\alias{fc2_2_select_delta_tau_BIC}
\alias{fc3_1_select_delta_tau_BIC}
\alias{fc3_2_select_delta_tau_BIC}
\alias{fc_1_select_delta_tau_CV}
\alias{fc_2_select_delta_tau_CV}
\alias{fr_select_delta_tau}
\alias{fr_partition_phi_N}
\alias{fr_onerep_tau}
\alias{fr_onerep_givenGroups_tau}
\alias{fr_onerep_givenThresholds_tau}
\alias{fc_ifCovers_vec}
\alias{fr_ifCover}
\alias{check_function}
\alias{VaRbacktestB}
\alias{fr_rolling_forecast}
\alias{fr2_rolling_forecast}
\alias{fr3_rolling_forecast_summary}
\title{Panel Quantile GARCH models under Homogeneity}
\description{
  Panel Quantile GARCH models under Homogeneity
}
\usage{
  fr_onerep_tau(tau, Y_N_T, if_group=c(TRUE, TRUE, TRUE), if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123_candid_lower=c(0,0,0), delta_123_candid_upper=c(2,2,2), testLength=NA);
  fr2_rolling_forecast(t0, tau, Y_N_T, if_group=c(TRUE, TRUE, TRUE), if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123_candid_lower=c(0,0,0), delta_123_candid_upper=c(2,2,2), testLength=NA, clsNum);
}
\arguments{
  \item{t0}{integer, moving window in rolling forecasting procedure}
  \item{tau}{decimal in (0,1), quantile level}
  \item{Y_N_T}{(N, T) matrix, panel data}
  \item{clsNum}{integer, number of core in parallel}
  \item{if_group}{logical vector, if there are homogeneous structures in the coefficient functions}
  \item{criterion}{character, default value is "BIC3", that is using the Bayesian information criterion (BIC) to choose the threshold for each group separately in binary segmentation}
  \item{delta_123_candid_lower}{vector, lower bound of grid search in selecting thresholds in binary segmentation}
  \item{delta_123_candid_upper}{vector, upper bound of grid search in selecting thresholds in binary segmentation}
  \item{if_mean}{logical, setting in selecting thresholds in binary segmentation}
  \item{tolerance}{decimal, setting in selecting thresholds in binary segmentation}
  \item{num_candid_onesplit}{integer, setting in selecting thresholds in binary segmentation}
  \item{testLength}{integer, valid only when criterion="CV"}
}
\value{
  \code{fr_onerep_tau} returns a vector, and \code{fr2_rolling_forecast} returns a list.
}
\details{
  With the preproccess in following "Examples", we can get the estimated coeﬀicients of the panel quantile GARCH model at 2.5\% quantile level using the three-stage estimation procedure.
}
\examples{
  library(PanelQGARCH)
  data("Y_N96_T3141")
  data("Initial_Estimation")
  
  tau <- 0.025 # quantile level
  N <- nrow(Y_N_T); T <- ncol(Y_N_T) # N and T
  if_group=c(TRUE, TRUE, TRUE) # the grouping structure varies across different parameters: \omega_i's, \alpha_'s and \beta_i's 
  criterion="BIC3" # use the Bayesian information criterion (BIC) to choose the threshold for each group separately in binary segmentation
  delta_123_candid_lower=c(0,0,0); delta_123_candid_upper=c(2,2,2) # selecting thresholds for binary segmentation via BIC-optimized grid search over [0, 2]
  if_mean=FALSE; tolerance=0.0001; num_candid_onesplit=6; testLength=NA # other setting in algorithm

  # est_tau2.5per_vec <- fr_onerep_tau(tau, Y_N_T, if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, testLength)
  
  ## Stage 1 - Preliminary estimation
  w_N_T <- fr_weights_NT(Y_N_T) # weights
  ### The initial estimation are obtained using R code by "Quantile autoregressive conditional heteroscedasticity"
  ### As this optimization algorithm has a certain degree of randomness, the results may vary slightly under different random seeds or software versions. To ensure consistency with the results reported in the paper, we provide the initial estimates in this package.
  ### phitilde_N_3 <- Initial_tau2.5per_N_3, that is the initial estimates in paper
  ### This package implements the functionality through the function fc_phitilde_N_tau().
  phitilde_N_3 <- fc_phitilde_N_tau(tau, Y_N_T, w_N_T) # \widetilde{\bm\phi}
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## Stage 2 - Homogeneity pursuit
  ### threshold selection
  w_N_T_train <- NA; phitilde_N_3_train <- NA
  if(criterion == "CV") { # initial estimates (training set)
    trainLength <- T - testLength; Y_N_T_train <- Y_N_T[,1:trainLength]; Y_N_T_test <- Y_N_T[,trainLength+c(1:testLength)]
    w_N_T_train <- fr_weights_NT(Y_N_T_train)
    phitilde_N_3_train <- fc_phitilde_N_tau(tau, Y_N_T_train, w_N_T_train)
  }
  deltahat_123 <- fr_select_delta_tau(if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, tau, Y_N_T, w_N_T, phitilde_N_3, testLength, w_N_T_train, phitilde_N_3_train) # (\widehat{\delta}_1, \widehat{\delta}_2, \widehat{\delta}_3)
  ### change points detection
  ghat_1 <- 1; changepoints_1 <- c(N); rank_1 <- rank(omegatilde_N)
  ghat_2 <- 1; changepoints_2 <- c(N); rank_2 <- rank(alphatilde_N)
  ghat_3 <- 1; changepoints_3 <- c(N); rank_3 <- rank(betatilde_N)
  if(if_group[1] == TRUE) {
    omegatilde_N_ascend <- sort(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
    delta_omega <- deltahat_123[1]; changepoints_omegatilde <- fc_changepoints(delta_omega, 1, N, omegatilde_N_ascend); num_changepoints_omegatilde <- length(changepoints_omegatilde) # (\widehat{h}_{1,1}, ..., \widehat{h}_{1,\widehat{g}_1})' and \widehat{g}_1 
    ghat_1 <- num_changepoints_omegatilde; changepoints_1 <- changepoints_omegatilde # \widehat{g}_1, and change points in {\widetilde{\omega}_{(i)}}_{i=1}^{N}
  }
  if(if_group[2] == TRUE) {
    alphatilde_N_ascend <- sort(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
    delta_alpha <- deltahat_123[2]; changepoints_alphatilde <- fc_changepoints(delta_alpha, 1, N, alphatilde_N_ascend); num_changepoints_alphatilde <- length(changepoints_alphatilde) # (\widehat{h}_{2,1}, ..., \widehat{h}_{2,\widehat{g}_1})' and \widehat{g}_2 
    ghat_2 <- num_changepoints_alphatilde; changepoints_2 <- changepoints_alphatilde # \widehat{g}_2, and change points in {\widetilde{\alpha}_{(i)}}_{i=1}^{N}
  }
  if(if_group[3] == TRUE) {
    betatilde_N_ascend <- sort(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
    delta_beta <- deltahat_123[3]; changepoints_betatilde <- fc_changepoints(delta_beta, 1, N, betatilde_N_ascend); num_changepoints_betatilde <- length(changepoints_betatilde) # (\widehat{h}_{3,1}, ..., \widehat{h}_{3,\widehat{g}_1})' and \widehat{g}_3 
    ghat_3 <- num_changepoints_betatilde; changepoints_3 <- changepoints_betatilde # \widehat{g}_3, and change points in {\widetilde{\beta}_{(i)}}_{i=1}^{N}
  }
  ### partition detection
  partition_phihat_N_3 <- fr_partition_phi_N(ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) # group labels of individuals wrt parameters
  rownames(partition_phihat_N_3) <- 1:N
  colnames(partition_phihat_N_3) <- c("label_omega", "label_alpha", "label_beta")
  ## Stage 3 - Final estimation
  thetahat <- fr_thetahat_tau(tau, Y_N_T, w_N_T, phitilde_N_3, ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) # \widehat{\bm\theta}
  phihat_N_3 <- fr_phi_N_homo(thetahat, ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) # \widehat{\bm\phi}
  rownames(phihat_N_3) <- 1:N
  colnames(phihat_N_3) <- c("omegahat", "alphahat", "betahat")
}
% \references{
  % See the documentation for Armadillo, and RcppArmadillo, for more details.
% }
\author{Wenyu Li}

