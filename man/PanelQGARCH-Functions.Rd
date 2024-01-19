\name{PanelQGARCH-Functions}
\alias{fc_weighs_NT}
\alias{fc_phitilde_N_tau}
\alias{fc_q_y_t0plus1}
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
\alias{fc_ifCovers_vec}
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
\alias{fr_weights_NT}
\alias{fr_select_delta_tau}
\alias{fr_partition_phi_N}
\alias{fr_onerep_tau}
\alias{fr_oracle_onerep_tau}
\alias{fr_onerep_givenGroups_tau}
\alias{fr_ifCover}
\alias{VaRbacktestB}
\alias{fr_rolling_forecast}
\alias{fr_rolling_forecast_summary_average}
\alias{fr2_rolling_forecast_summary_average}
\title{Panel Quantile GARCH models under Homogeneity}
\description{
  Panel Quantile GARCH models under Homogeneity
}
\usage{
  fr_onerep_tau(tau, Y_N_T);
  fr_rolling_forecast(t0, tau, Y_N_T, clsNum);
}
\arguments{
  \item{t0}{integer, moving window in rolling forecasting procedure}
  \item{tau}{decimal in (0,1), quantile level}
  \item{Y_N_T}{(N, T) matrix, panel data}
  \item{clsNum}{integer, number of core in parallel}
}
\value{
  \code{fr_onerep_tau} returns a vector, and \code{fr_rolling_forecast} returns a matrix.
}
\details{
  With the preproccess in following "Examples", we can get the estimated coeﬀicients of the panel quantile GARCH model at 5\% quantile level using the three-stage estimation procedure.
}
\examples{
  library(PanelQGARCH)
  data("Y_N100_T1761")
  est_tau5per_vec <- fr_onerep_tau(tau=0.05, Y_N_T)
}
% \references{
  % See the documentation for Armadillo, and RcppArmadillo, for more details.
% }
\author{Wenyu Li}

