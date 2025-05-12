\name{PanelQGARCHwithFactor-Functions}
\alias{check_function}
\alias{psi_function}
\alias{h_compute}
\alias{arch_fft}
\alias{Init_par}
\alias{Loss_QR}
\alias{Loss_QR_gr}
\alias{fit1_optim}
\alias{fr_phitilde_i_tau}
\alias{fit1_optim_grid}
\alias{fc_weighs_NT}
\alias{fc_phitilde_N_tau}
\alias{fr_weights_NT}
\alias{fr_phi_l_N_homo}
\alias{fr_loss_panelQR}
\alias{fr2_loss_panelQR}
\alias{fr_theta_l_initial_IN_optim}
\alias{fr_thetahat_tau}
\alias{fr_thetahat_l_tau}
\alias{fr3_BIC_delta_l_tau}
\alias{fc_rank_mat}
\alias{fc_phi_N_homo}
\alias{fc_calculate_segment_means}
\alias{fc_theta_initial_IN_optim}
\alias{fc_Delta_ijh}
\alias{fc_changepoints}
\alias{fc_g_changepoints_vec}
\alias{fc_loss_panelQR}
\alias{fc3_1_select_delta_tau_BIC}
\alias{fc3_2_select_delta_tau_BIC}
\alias{fr_select_delta_tau}
\alias{fr_partition_phi_l_N}
\alias{fc_partition_phi_N}
\alias{fr_onerep_tau}
\title{Panel Quantile GARCH models with factors under Homogeneity}
\description{
  Panel Quantile GARCH models with factors under Homogeneity
}
\usage{
  fr_onerep_tau(tau, Y_N_T, Factors_r_T, if_group, if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123r_candid_lower, delta_123r_candid_upper);
}
\arguments{
  \item{tau}{decimal in (0,1), quantile level}
  \item{Y_N_T}{(N, T) matrix, panel data}
  \item{Factors_r_T}{(r, T) matrix, factors data}
  \item{if_group}{logical vector, if there are homogeneous structures in the coefficient functions}
  \item{criterion}{character, default value is "BIC3", that is using the Bayesian information criterion (BIC) to choose the threshold for each group separately in binary segmentation}
  \item{delta_123r_candid_lower}{vector, lower bound of grid search in selecting thresholds in binary segmentation}
  \item{delta_123r_candid_upper}{vector, upper bound of grid search in selecting thresholds in binary segmentation}
  \item{if_mean}{logical, setting in selecting thresholds in binary segmentation}
  \item{tolerance}{decimal, setting in selecting thresholds in binary segmentation}
  \item{num_candid_onesplit}{integer, setting in selecting thresholds in binary segmentation}
}
\value{
  \code{fr_onerep_tau} returns a vector.
}
\details{
  With the preproccess in following "Examples", we can get the estimated coeﬀicients of the panel quantile GARCH model with factors at 2.5\% quantile level using the three-stage estimation procedure.
}
\examples{
  library(PanelQGARCHwithFactor)
  data("Y_N30_T3142")
  data("Factors_r5_T3142")
  data("Initial_Estimation_factor")
  
  tau <- 0.025 # quantile level
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  r <- nrow(Factors_r_T)
  
  if_group=rep(TRUE,3+r) # the grouping structure varies across different parameters: \omega_i's, \alpha_'s and \beta_i's 
  criterion="BIC3" # use the Bayesian information criterion (BIC) to choose the threshold for each group separately in binary segmentation
  delta_123r_candid_lower=rep(0,3+r); delta_123r_candid_upper=rep(2,3+r) # selecting thresholds for binary segmentation via BIC-optimized grid search over [0, 2]
  if_mean=FALSE; tolerance=0.0001; num_candid_onesplit=6; testLength=NA # other setting in algorithm

  # est_tau2.5per_vec <- fr_onerep_tau(tau, Y_N_T, Factors_r_T, if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123r_candid_lower, delta_123r_candid_upper)
  
  ## Stage 1 - Initial estimates
  w_N_T <- fr_weights_NT(Y_N_T) # weights
  ### The initial estimation are obtained by extending R code in "Quantile autoregressive conditional heteroscedasticity"
  ### As this optimization algorithm has a certain degree of randomness, the results may vary slightly under different random seeds or software versions. To ensure consistency with the results reported in the paper, we provide the initial estimates in this package.
  ### phitilde_N_3plusr <- Initial_tau2.5per_N_8, that is the initial estimates in paper
  ### This package implements the functionality through the function fc_phitilde_N_tau().
  phitilde_N_3plusr <- fc_phitilde_N_tau(tau, Y_N_T, Factors_r_T, w_N_T) # \widetilde{\bm\phi}
  ## Stage 2 - Change points detection by binary segmentation algorithm
  ### threshold selection
  rank_mat <- fc_rank_mat(phitilde_N_3plusr)
  deltahat_123r <- fr_select_delta_tau(if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123r_candid_lower, delta_123r_candid_upper, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat) # {\widehat{\delta}}_{\ell=1}^{8}
  ### change points detection
  ghat_changepoints_vec <- fc_g_changepoints_vec(phitilde_N_3plusr, deltahat_123r)
  ghat_vec <- ghat_changepoints_vec[1:(3+r)] # {\widehat{g}_\ell}_{\ell=1}^{8}
  changepoints_vec <- ghat_changepoints_vec[-c(1:(3+r))] # change points
  ### partition detection
  partition_phihat_N_3plusr <- fc_partition_phi_N(ghat_vec, changepoints_vec, rank_mat) # group labels of individuals wrt parameters
  rownames(partition_phihat_N_3plusr) <- 1:N
  colnames(partition_phihat_N_3plusr) <- c("label_omega", "label_alpha", "label_beta", paste0("label_lambda", 1:5))
  ## Stage 3 - Final estimates
  thetahat <- fr_thetahat_tau(tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, ghat_vec, changepoints_vec, rank_mat) # \widehat{\bm\theta}
  phihat_N_3plusr <- fc_phi_N_homo(thetahat, ghat_vec, changepoints_vec, rank_mat) # \widehat{\bm\phi}
  rownames(phihat_N_3plusr) <- 1:N
  colnames(phihat_N_3plusr) <- c("omegahat", "alphahat", "betahat",  paste0("lambda", 1:5, "hat"))
}
% \references{
  % See the documentation for Armadillo, and RcppArmadillo, for more details.
% }
\author{Wenyu Li}

