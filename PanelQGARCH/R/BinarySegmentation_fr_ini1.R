### selection of thresholds \delta_1, \delta_2, \delta_3
fr_phi_N_homo <- function(theta_vec, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) {
  # transform the parameter vector \bm\theta = (\omega_{(1)}, ..., \omega_{(g_1)}, \alpha_{(1)}, ..., \alpha_{(g_2)}, \beta_{(1)}, ..., \beta_{(g_3)})' 
  # to the N \times 3 parameter matrix (\bm{\phi}_1, ..., \bm{\phi}_N)' under the homogeneity structure
  omega_g1_ascend <- theta_vec[1:g_1] # (\omega_{(1)}, ..., \omega_{(g_1)})' 
  alpha_g2_ascend <- theta_vec[g_1+(1:g_2)] # (\alpha_{(1)}, ..., \alpha_{(g_2)})'
  beta_g3_ascend <- theta_vec[g_1+g_2+(1:g_3)] # (\beta_{(1)}, ..., \beta_{(g_3)})'
  omega_N_ascend <- rep(omega_g1_ascend, diff(c(0, changepoints_1))) # (\omega_{(1)}, ..., \omega_{(N)})' 
  alpha_N_ascend <- rep(alpha_g2_ascend, diff(c(0, changepoints_2))) # (\alpha_{(1)}, ..., \alpha_{(N)})' 
  beta_N_ascend <- rep(beta_g3_ascend, diff(c(0, changepoints_3))) # (\beta_{(1)}, ..., \beta_{(N)})' 
  omega_N <- omega_N_ascend[rank_1]
  alpha_N <- alpha_N_ascend[rank_2]
  beta_N <- beta_N_ascend[rank_3]
  phi_N_3 <- cbind(omega_N, alpha_N, beta_N)
  return(phi_N_3)
}
# #### example
# omegatilde_N <- c(5,5,6,6,6,8,8,2,9,9,9,9,9,3,3,3)
# theta <- c(2,3,5,6,8,9)
# ##### test
# ###### detection
# omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
# delta_omega <- 0.8; changepoints_omegatilde <- fc_changepoints(delta_omega, 1, length(omegatilde_N), omegatilde_N_ascend); num_changepoints_omegatilde <- length(changepoints_omegatilde) # (\widehat{h}_{1,1}, ..., \widehat{h}_{1,\widehat{g}_1})' and \widehat{g}_1
# omegatilde_N_ascend; changepoints_omegatilde; num_changepoints_omegatilde
# ###### transformation
# g_1 <- num_changepoints_omegatilde; changepoints_1 <- changepoints_omegatilde; rank_1 <- rank_omegatilde_N
# omega_g1_ascend <- theta[1:g_1]
# omega_N_ascend <- rep(omega_g1_ascend, diff(c(0, changepoints_1))) # (\omega_{(1)}, ..., \omega_{(N)})'
# omega_N <- omega_N_ascend[rank_1]
# omega_N_ascend; omegatilde_N_ascend
# omega_N; omegatilde_N

fr_loss_panelQR <- function(tau, theta_vec, Y_N_T, w_N_T, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) {
  # the loss function of the final estimator 
  phi_N_3 <- fr_phi_N_homo(theta_vec, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  loss_PanelQR <- fc_loss_panelQR(tau, phi_N_3, Y_N_T, w_N_T)
  return(loss_PanelQR)
}

fr2_loss_panelQR <- function(tau, theta_vec_omega, theta_vec_alpha, theta_vec_beta, Y_N_T, w_N_T, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) {
  # the loss function of the final estimator 
  theta_vec <- c(theta_vec_omega, theta_vec_alpha, theta_vec_beta)
  phi_N_3 <- fr_phi_N_homo(theta_vec, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  loss_PanelQR <- fc_loss_panelQR(tau, phi_N_3, Y_N_T, w_N_T)
  return(loss_PanelQR)
}

fr_theta_initial_IN_optim <- function(phiini_N_3, changepoints_1, changepoints_2, changepoints_3) {
  # initial value in optim() for estimating \widehat\bm\theta
  omegaini_N <- phiini_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphaini_N <- phiini_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betaini_N <- phiini_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  omegaini_N_ascend <- sort(omegaini_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphaini_N_ascend <- sort(alphaini_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betaini_N_ascend <- sort(betaini_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  theta_ini <- c(omegaini_N_ascend[changepoints_1], alphaini_N_ascend[changepoints_2], betaini_N_ascend[changepoints_3])
  return(theta_ini)
}
# fr_theta_initial_IN_optim <- function(phiini_N_3, changepoints_1, changepoints_2, changepoints_3) {
#   # initial value in optim() for estimating \widehat\bm\theta
#   omegaini_N <- phiini_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
#   alphaini_N <- phiini_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
#   betaini_N <- phiini_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
#   omegaini_N_ascend <- sort(omegaini_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
#   alphaini_N_ascend <- sort(alphaini_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
#   betaini_N_ascend <- sort(betaini_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
#   # theta_ini <- c(omegaini_N_ascend[changepoints_1], alphaini_N_ascend[changepoints_2], betaini_N_ascend[changepoints_3])
#   omegaINtheta_ini <- fc_calculate_segment_means(omegaini_N_ascend, c(0, changepoints_1))
#   alphaINtheta_ini <- fc_calculate_segment_means(alphaini_N_ascend, c(0, changepoints_2))
#   betaINtheta_ini <- fc_calculate_segment_means(betaini_N_ascend, c(0, changepoints_3))
#   theta_ini <- c(omegaINtheta_ini, alphaINtheta_ini, betaINtheta_ini)
#   return(theta_ini)
# }

fr_thetahat_tau <- function(tau, Y_N_T, w_N_T, phitilde_N_3, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) {
  # the final estimator \widehat\bm\theta = (\widehat\omega_{(1)}, ..., \widehat\omega_{(\widehat{g}_1)}, \widehat\alpha_{(1)}, ..., \widehat\alpha_{(\widehat{g}_2)}, \widehat\beta_{(1)}, ..., \widehat\beta_{(\widehat{g}_3)})', which is obtained using the threshold \delta
  theta_ini <- fr_theta_initial_IN_optim(phitilde_N_3, changepoints_1, changepoints_2, changepoints_3)
  fit_thetahat <- optim(par=theta_ini, fn=fr_loss_panelQR, method="L-BFGS-B", lower=c(rep(-Inf,g_1+g_2),rep(0,g_3)), upper=c(rep(Inf,g_1+g_2),rep(1-1e-3,g_3)), tau=tau, Y_N_T=Y_N_T, w_N_T=w_N_T, g_1=g_1, g_2=g_2, g_3=g_3, changepoints_1=changepoints_1, changepoints_2=changepoints_2, changepoints_3=changepoints_3, rank_1=rank_1, rank_2=rank_2, rank_3=rank_3) 
  thetahat <- fit_thetahat$par
  return(thetahat)
}

fr_thetahat_omega_tau <- function(tau, Y_N_T, w_N_T, phitilde_N_3, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) {
  # (\widehat\omega_{(1)}, ..., \widehat\omega_{(\widehat{g}_1)})' when given (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)}, \widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  theta_ini <- fr_theta_initial_IN_optim(phitilde_N_3, changepoints_1, changepoints_2, changepoints_3)
  fit_thetahat_omega <- optim(par=theta_ini[1:g_1], fn=fr2_loss_panelQR, method="L-BFGS-B", lower=rep(-Inf,g_1), upper=rep(Inf,g_1), tau=tau, theta_vec_alpha=theta_ini[g_1+(1:g_2)], theta_vec_beta=theta_ini[g_1+g_2+(1:g_3)], Y_N_T=Y_N_T, w_N_T=w_N_T, g_1=g_1, g_2=g_2, g_3=g_3, changepoints_1=changepoints_1, changepoints_2=changepoints_2, changepoints_3=changepoints_3, rank_1=rank_1, rank_2=rank_2, rank_3=rank_3) 
  thetahat_omega <- fit_thetahat_omega$par
  return(thetahat_omega)
}

fr_thetahat_alpha_tau <- function(tau, Y_N_T, w_N_T, phitilde_N_3, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) {
  # (\widehat\alpha_{(1)}, ..., \widehat\alpha_{(\widehat{g}_2)})' when given (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)}, \widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  theta_ini <- fr_theta_initial_IN_optim(phitilde_N_3, changepoints_1, changepoints_2, changepoints_3)
  fit_thetahat_alpha <- optim(par=theta_ini[g_1+(1:g_2)], fn=fr2_loss_panelQR, method="L-BFGS-B", lower=rep(-Inf,g_2), upper=rep(Inf,g_2), tau=tau, theta_vec_omega=theta_ini[1:g_1], theta_vec_beta=theta_ini[g_1+g_2+(1:g_3)], Y_N_T=Y_N_T, w_N_T=w_N_T, g_1=g_1, g_2=g_2, g_3=g_3, changepoints_1=changepoints_1, changepoints_2=changepoints_2, changepoints_3=changepoints_3, rank_1=rank_1, rank_2=rank_2, rank_3=rank_3) 
  thetahat_alpha <- fit_thetahat_alpha$par
  return(thetahat_alpha)
}

fr_thetahat_beta_tau <- function(tau, Y_N_T, w_N_T, phitilde_N_3, g_1, g_2, g_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3) {
  # (\widehat\beta_{(1)}, ..., \widehat\beta_{(\widehat{g}_3)})' when given (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)}, \widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  theta_ini <- fr_theta_initial_IN_optim(phitilde_N_3, changepoints_1, changepoints_2, changepoints_3)
  fit_thetahat_beta <- optim(par=theta_ini[g_1+g_2+(1:g_3)], fn=fr2_loss_panelQR, method="L-BFGS-B", lower=rep(0,g_3), upper=rep(1-1e-3,g_3), tau=tau, theta_vec_omega=theta_ini[1:g_1], theta_vec_alpha=theta_ini[g_1+(1:g_2)], Y_N_T=Y_N_T, w_N_T=w_N_T, g_1=g_1, g_2=g_2, g_3=g_3, changepoints_1=changepoints_1, changepoints_2=changepoints_2, changepoints_3=changepoints_3, rank_1=rank_1, rank_2=rank_2, rank_3=rank_3) 
  thetahat_beta <- fit_thetahat_beta$par
  return(thetahat_beta)
}

fr_BIC_delta_tau <- function(delta_123, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_omega <- delta_123[1]; changepoints_omegatilde <- fc_changepoints(delta_omega, 1, N, omegatilde_N_ascend); num_changepoints_omegatilde <- length(changepoints_omegatilde) # (\widehat{h}_{1,1}, ..., \widehat{h}_{1,\widehat{g}_1})' and \widehat{g}_1 
  delta_alpha <- delta_123[2]; changepoints_alphatilde <- fc_changepoints(delta_alpha, 1, N, alphatilde_N_ascend); num_changepoints_alphatilde <- length(changepoints_alphatilde) # (\widehat{h}_{2,1}, ..., \widehat{h}_{2,\widehat{g}_1})' and \widehat{g}_2 
  delta_beta <- delta_123[3]; changepoints_betatilde <- fc_changepoints(delta_beta, 1, N, betatilde_N_ascend); num_changepoints_betatilde <- length(changepoints_betatilde) # (\widehat{h}_{3,1}, ..., \widehat{h}_{3,\widehat{g}_1})' and \widehat{g}_3 
  ghat_1 <- num_changepoints_omegatilde; changepoints_1 <- changepoints_omegatilde; rank_1 <- rank_omegatilde_N
  ghat_2 <- num_changepoints_alphatilde; changepoints_2 <- changepoints_alphatilde; rank_2 <- rank_alphatilde_N
  ghat_3 <- num_changepoints_betatilde; changepoints_3 <- changepoints_betatilde; rank_3 <- rank_betatilde_N
  ## final estimates
  thetahat <- fr_thetahat_tau(tau, Y_N_T, w_N_T, phitilde_N_3, ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  # BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  BIC <- 2*N*log(loss_panelQR/N) + length(thetahat)*log(N)
  return(BIC)
}

fr_CV_delta_tau <- function(delta_123, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train) {
  # determine the threshold \bm\delta = (\delta_1, \delta_2, \delta_3)' via cross-validation
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  trainLength <- T - testLength; Y_N_T_train <- Y_N_T[,1:trainLength]; Y_N_T_test <- Y_N_T[,trainLength+c(1:testLength)]
  ## initial estimates (training set)
  # w_N_T_train <- fr_weights_NT(Y_N_T_train)
  # phitilde_N_3_train <- fc_phitilde_N_tau(tau, Y_N_T_train, w_N_T_train)
  omegatilde_N <- phitilde_N_3_train[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3_train[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3_train[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection (training set)
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_omega <- delta_123[1]; changepoints_omegatilde <- fc_changepoints(delta_omega, 1, N, omegatilde_N_ascend); num_changepoints_omegatilde <- length(changepoints_omegatilde) # (\widehat{h}_{1,1}, ..., \widehat{h}_{1,\widehat{g}_1})' and \widehat{g}_1 
  delta_alpha <- delta_123[2]; changepoints_alphatilde <- fc_changepoints(delta_alpha, 1, N, alphatilde_N_ascend); num_changepoints_alphatilde <- length(changepoints_alphatilde) # (\widehat{h}_{2,1}, ..., \widehat{h}_{2,\widehat{g}_1})' and \widehat{g}_2 
  delta_beta <- delta_123[3]; changepoints_betatilde <- fc_changepoints(delta_beta, 1, N, betatilde_N_ascend); num_changepoints_betatilde <- length(changepoints_betatilde) # (\widehat{h}_{3,1}, ..., \widehat{h}_{3,\widehat{g}_1})' and \widehat{g}_3 
  ghat_1 <- num_changepoints_omegatilde; changepoints_1 <- changepoints_omegatilde; rank_1 <- rank_omegatilde_N
  ghat_2 <- num_changepoints_alphatilde; changepoints_2 <- changepoints_alphatilde; rank_2 <- rank_alphatilde_N
  ghat_3 <- num_changepoints_betatilde; changepoints_3 <- changepoints_betatilde; rank_3 <- rank_betatilde_N
  ## final estimates (training set)
  thetahat <- fr_thetahat_tau(tau, Y_N_T_train, w_N_T_train, phitilde_N_3_train, ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  ## CV (validation set)
  w_N_T_test <- fr_weights_NT(Y_N_T_test)
  prediction_error <- fr_loss_panelQR(tau, thetahat, Y_N_T_test, w_N_T_test, ghat_1, ghat_2, ghat_3, changepoints_1, changepoints_2, changepoints_3, rank_1, rank_2, rank_3)
  return(prediction_error)
}

fr1_BIC_delta_1_tau <- function(delta_1, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \delta_1 using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_omega <- delta_1; changepoints_omegatilde <- fc_changepoints(delta_omega, 1, N, omegatilde_N_ascend); num_changepoints_omegatilde <- length(changepoints_omegatilde) # (\widehat{h}_{1,1}, ..., \widehat{h}_{1,\widehat{g}_1})' and \widehat{g}_1 
  ghat_1 <- num_changepoints_omegatilde; changepoints_1 <- changepoints_omegatilde; rank_1 <- rank_omegatilde_N
  ## final estimates
  thetahat_omega <- fr_thetahat_omega_tau(tau, Y_N_T, w_N_T, phitilde_N_3, ghat_1, N, N, changepoints_1, c(1:N), c(1:N), rank_1, rank_alphatilde_N, rank_betatilde_N)
  thetahat <- c(thetahat_omega, alphatilde_N_ascend, betatilde_N_ascend)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, ghat_1, N, N, changepoints_1, c(1:N), c(1:N), rank_1, rank_alphatilde_N, rank_betatilde_N)
  BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  return(BIC)
}

fr1_BIC_delta_2_tau <- function(delta_2, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \delta_2 using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_alpha <- delta_2; changepoints_alphatilde <- fc_changepoints(delta_alpha, 1, N, alphatilde_N_ascend); num_changepoints_alphatilde <- length(changepoints_alphatilde) # (\widehat{h}_{2,1}, ..., \widehat{h}_{2,\widehat{g}_1})' and \widehat{g}_2 
  ghat_2 <- num_changepoints_alphatilde; changepoints_2 <- changepoints_alphatilde; rank_2 <- rank_alphatilde_N
  ## final estimates
  thetahat_alpha <- fr_thetahat_alpha_tau(tau, Y_N_T, w_N_T, phitilde_N_3, N, ghat_2, N, c(1:N), changepoints_2, c(1:N), rank_omegatilde_N, rank_2, rank_betatilde_N)
  thetahat <- c(omegatilde_N_ascend, thetahat_alpha, betatilde_N_ascend)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, N, ghat_2, N, c(1:N), changepoints_2, c(1:N), rank_omegatilde_N, rank_2, rank_betatilde_N)
  BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  return(BIC)
}

fr1_BIC_delta_3_tau <- function(delta_3, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \delta_3 using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_beta <- delta_3; changepoints_betatilde <- fc_changepoints(delta_beta, 1, N, betatilde_N_ascend); num_changepoints_betatilde <- length(changepoints_betatilde) # (\widehat{h}_{3,1}, ..., \widehat{h}_{3,\widehat{g}_1})' and \widehat{g}_3 
  ghat_3 <- num_changepoints_betatilde; changepoints_3 <- changepoints_betatilde; rank_3 <- rank_betatilde_N
  ## final estimates
  thetahat_beta <- fr_thetahat_beta_tau(tau, Y_N_T, w_N_T, phitilde_N_3, N, N, ghat_3, c(1:N), c(1:N), changepoints_3, rank_omegatilde_N, rank_alphatilde_N, rank_3)
  thetahat <- c(omegatilde_N_ascend, alphatilde_N_ascend, thetahat_beta)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, N, N, ghat_3, c(1:N), c(1:N), changepoints_3, rank_omegatilde_N, rank_alphatilde_N, rank_3)
  BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  return(BIC)
}

fr2_BIC_delta_1_tau <- function(delta_1, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \delta_1 using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_omega <- delta_1; changepoints_omegatilde <- fc_changepoints(delta_omega, 1, N, omegatilde_N_ascend); num_changepoints_omegatilde <- length(changepoints_omegatilde) # (\widehat{h}_{1,1}, ..., \widehat{h}_{1,\widehat{g}_1})' and \widehat{g}_1 
  ghat_1 <- num_changepoints_omegatilde; changepoints_1 <- changepoints_omegatilde; rank_1 <- rank_omegatilde_N
  ## final estimates
  thetahat_omega <- fr_thetahat_omega_tau(tau, Y_N_T, w_N_T, phitilde_N_3, ghat_1, N, N, changepoints_1, c(1:N), c(1:N), rank_1, rank_alphatilde_N, rank_betatilde_N)
  thetahat <- c(thetahat_omega, alphatilde_N_ascend, betatilde_N_ascend)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, ghat_1, N, N, changepoints_1, c(1:N), c(1:N), rank_1, rank_alphatilde_N, rank_betatilde_N)
  # BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  BIC <- 2*N*log(loss_panelQR/N) + ghat_1*log(N)
  return(BIC)
}

fr2_BIC_delta_2_tau <- function(delta_2, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \delta_2 using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_alpha <- delta_2; changepoints_alphatilde <- fc_changepoints(delta_alpha, 1, N, alphatilde_N_ascend); num_changepoints_alphatilde <- length(changepoints_alphatilde) # (\widehat{h}_{2,1}, ..., \widehat{h}_{2,\widehat{g}_1})' and \widehat{g}_2 
  ghat_2 <- num_changepoints_alphatilde; changepoints_2 <- changepoints_alphatilde; rank_2 <- rank_alphatilde_N
  ## final estimates
  thetahat_alpha <- fr_thetahat_alpha_tau(tau, Y_N_T, w_N_T, phitilde_N_3, N, ghat_2, N, c(1:N), changepoints_2, c(1:N), rank_omegatilde_N, rank_2, rank_betatilde_N)
  thetahat <- c(omegatilde_N_ascend, thetahat_alpha, betatilde_N_ascend)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, N, ghat_2, N, c(1:N), changepoints_2, c(1:N), rank_omegatilde_N, rank_2, rank_betatilde_N)
  # BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  BIC <- 2*N*log(loss_panelQR/N) + ghat_2*log(N)
  return(BIC)
}

fr2_BIC_delta_3_tau <- function(delta_3, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \delta_3 using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_beta <- delta_3; changepoints_betatilde <- fc_changepoints(delta_beta, 1, N, betatilde_N_ascend); num_changepoints_betatilde <- length(changepoints_betatilde) # (\widehat{h}_{3,1}, ..., \widehat{h}_{3,\widehat{g}_1})' and \widehat{g}_3 
  ghat_3 <- num_changepoints_betatilde; changepoints_3 <- changepoints_betatilde; rank_3 <- rank_betatilde_N
  ## final estimates
  thetahat_beta <- fr_thetahat_beta_tau(tau, Y_N_T, w_N_T, phitilde_N_3, N, N, ghat_3, c(1:N), c(1:N), changepoints_3, rank_omegatilde_N, rank_alphatilde_N, rank_3)
  thetahat <- c(omegatilde_N_ascend, alphatilde_N_ascend, thetahat_beta)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, N, N, ghat_3, c(1:N), c(1:N), changepoints_3, rank_omegatilde_N, rank_alphatilde_N, rank_3)
  # BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  BIC <- 2*N*log(loss_panelQR/N) + ghat_3*log(N)
  return(BIC)
}

fr3_BIC_delta_1_tau <- function(delta_1, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \delta_1 using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_omega <- delta_1; changepoints_omegatilde <- fc_changepoints(delta_omega, 1, N, omegatilde_N_ascend); num_changepoints_omegatilde <- length(changepoints_omegatilde) # (\widehat{h}_{1,1}, ..., \widehat{h}_{1,\widehat{g}_1})' and \widehat{g}_1 
  ghat_1 <- num_changepoints_omegatilde; changepoints_1 <- changepoints_omegatilde; rank_1 <- rank_omegatilde_N
  ## final estimates
  thetahat_omega <- fr_thetahat_omega_tau(tau, Y_N_T, w_N_T, phitilde_N_3, ghat_1, N, N, changepoints_1, c(1:N), c(1:N), rank_1, rank_alphatilde_N, rank_betatilde_N)
  thetahat <- c(thetahat_omega, alphatilde_N_ascend, betatilde_N_ascend)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, ghat_1, N, N, changepoints_1, c(1:N), c(1:N), rank_1, rank_alphatilde_N, rank_betatilde_N)
  # BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  # BIC <- 2*N*log(loss_panelQR/N) + ghat_1*log(N)
  BIC <- loss_panelQR + ghat_1*log(N)
  return(BIC)
}

fr3_BIC_delta_2_tau <- function(delta_2, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \delta_2 using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_alpha <- delta_2; changepoints_alphatilde <- fc_changepoints(delta_alpha, 1, N, alphatilde_N_ascend); num_changepoints_alphatilde <- length(changepoints_alphatilde) # (\widehat{h}_{2,1}, ..., \widehat{h}_{2,\widehat{g}_1})' and \widehat{g}_2 
  ghat_2 <- num_changepoints_alphatilde; changepoints_2 <- changepoints_alphatilde; rank_2 <- rank_alphatilde_N
  ## final estimates
  thetahat_alpha <- fr_thetahat_alpha_tau(tau, Y_N_T, w_N_T, phitilde_N_3, N, ghat_2, N, c(1:N), changepoints_2, c(1:N), rank_omegatilde_N, rank_2, rank_betatilde_N)
  thetahat <- c(omegatilde_N_ascend, thetahat_alpha, betatilde_N_ascend)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, N, ghat_2, N, c(1:N), changepoints_2, c(1:N), rank_omegatilde_N, rank_2, rank_betatilde_N)
  # BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  # BIC <- 2*N*log(loss_panelQR/N) + ghat_2*log(N)
  BIC <- loss_panelQR + ghat_2*log(N)
  return(BIC)
}

fr3_BIC_delta_3_tau <- function(delta_3, tau, Y_N_T, w_N_T, phitilde_N_3) {
  # determine the threshold \delta_3 using BIC
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  ## initial estimates
  omegatilde_N <- phitilde_N_3[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_beta <- delta_3; changepoints_betatilde <- fc_changepoints(delta_beta, 1, N, betatilde_N_ascend); num_changepoints_betatilde <- length(changepoints_betatilde) # (\widehat{h}_{3,1}, ..., \widehat{h}_{3,\widehat{g}_1})' and \widehat{g}_3 
  ghat_3 <- num_changepoints_betatilde; changepoints_3 <- changepoints_betatilde; rank_3 <- rank_betatilde_N
  ## final estimates
  thetahat_beta <- fr_thetahat_beta_tau(tau, Y_N_T, w_N_T, phitilde_N_3, N, N, ghat_3, c(1:N), c(1:N), changepoints_3, rank_omegatilde_N, rank_alphatilde_N, rank_3)
  thetahat <- c(omegatilde_N_ascend, alphatilde_N_ascend, thetahat_beta)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, w_N_T, N, N, ghat_3, c(1:N), c(1:N), changepoints_3, rank_omegatilde_N, rank_alphatilde_N, rank_3)
  # BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  # BIC <- 2*N*log(loss_panelQR/N) + ghat_3*log(N)
  BIC <- loss_panelQR + ghat_3*log(N)
  return(BIC)
}

fr_CV_delta_1_tau <- function(delta_1, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train) {
  # determine the threshold \delta_1 via cross-validation
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  trainLength <- T - testLength; Y_N_T_train <- Y_N_T[,1:trainLength]; Y_N_T_test <- Y_N_T[,trainLength+c(1:testLength)]
  ## initial estimates (training set)
  # w_N_T_train <- fr_weights_NT(Y_N_T_train)
  # phitilde_N_3_train <- fc_phitilde_N_tau(tau, Y_N_T_train, w_N_T_train)
  omegatilde_N <- phitilde_N_3_train[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3_train[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3_train[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection (training set)
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_omega <- delta_1; changepoints_omegatilde <- fc_changepoints(delta_omega, 1, N, omegatilde_N_ascend); num_changepoints_omegatilde <- length(changepoints_omegatilde) # (\widehat{h}_{1,1}, ..., \widehat{h}_{1,\widehat{g}_1})' and \widehat{g}_1 
  ghat_1 <- num_changepoints_omegatilde; changepoints_1 <- changepoints_omegatilde; rank_1 <- rank_omegatilde_N
  ## final estimates (training set)
  thetahat_omega <- fr_thetahat_omega_tau(tau, Y_N_T_train, w_N_T_train, phitilde_N_3_train, ghat_1, N, N, changepoints_1, c(1:N), c(1:N), rank_1, rank_alphatilde_N, rank_betatilde_N)
  thetahat <- c(thetahat_omega, alphatilde_N_ascend, betatilde_N_ascend)
  ## CV (validation set)
  w_N_T_test <- fr_weights_NT(Y_N_T_test)
  prediction_error <- fr_loss_panelQR(tau, thetahat, Y_N_T_test, w_N_T_test, ghat_1, N, N, changepoints_1, c(1:N), c(1:N), rank_1, rank_alphatilde_N, rank_betatilde_N)
  return(prediction_error)
}

fr_CV_delta_2_tau <- function(delta_2, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train) {
  # determine the threshold \delta_2 via cross-validation
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  trainLength <- T - testLength; Y_N_T_train <- Y_N_T[,1:trainLength]; Y_N_T_test <- Y_N_T[,trainLength+c(1:testLength)]
  ## initial estimates (training set)
  # w_N_T_train <- fr_weights_NT(Y_N_T_train)
  # phitilde_N_3_train <- fc_phitilde_N_tau(tau, Y_N_T_train, w_N_T_train)
  omegatilde_N <- phitilde_N_3_train[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3_train[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3_train[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection (training set)
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_alpha <- delta_2; changepoints_alphatilde <- fc_changepoints(delta_alpha, 1, N, alphatilde_N_ascend); num_changepoints_alphatilde <- length(changepoints_alphatilde) # (\widehat{h}_{2,1}, ..., \widehat{h}_{2,\widehat{g}_1})' and \widehat{g}_2 
  ghat_2 <- num_changepoints_alphatilde; changepoints_2 <- changepoints_alphatilde; rank_2 <- rank_alphatilde_N
  ## final estimates (training set)
  thetahat_alpha <- fr_thetahat_alpha_tau(tau, Y_N_T_train, w_N_T_train, phitilde_N_3_train, N, ghat_2, N, c(1:N), changepoints_2, c(1:N), rank_omegatilde_N, rank_2, rank_betatilde_N)
  thetahat <- c(omegatilde_N_ascend, thetahat_alpha, betatilde_N_ascend)
  ## CV (validation set)
  w_N_T_test <- fr_weights_NT(Y_N_T_test)
  prediction_error <- fr_loss_panelQR(tau, thetahat, Y_N_T_test, w_N_T_test, N, ghat_2, N, c(1:N), changepoints_2, c(1:N), rank_omegatilde_N, rank_2, rank_betatilde_N)
  return(prediction_error)
}

fr_CV_delta_3_tau <- function(delta_3, tau, Y_N_T, testLength, w_N_T_train, phitilde_N_3_train) {
  # determine the threshold \delta_3 via cross-validation
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  trainLength <- T - testLength; Y_N_T_train <- Y_N_T[,1:trainLength]; Y_N_T_test <- Y_N_T[,trainLength+c(1:testLength)]
  ## initial estimates (training set)
  # w_N_T_train <- fr_weights_NT(Y_N_T_train)
  # phitilde_N_3_train <- fc_phitilde_N_tau(tau, Y_N_T_train, w_N_T_train)
  omegatilde_N <- phitilde_N_3_train[,1] # (\widetilde\omega_1, ..., \widetilde\omega_N)'
  alphatilde_N <- phitilde_N_3_train[,2] # (\widetilde\alpha_1, ..., \widetilde\alpha_N)'
  betatilde_N <- phitilde_N_3_train[,3] # (\widetilde\beta_1, ..., \widetilde\beta_N)'
  ## change points detection (training set)
  omegatilde_N_ascend <- sort(omegatilde_N); rank_omegatilde_N <- rank(omegatilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})'
  alphatilde_N_ascend <- sort(alphatilde_N); rank_alphatilde_N <- rank(alphatilde_N) # (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})'
  betatilde_N_ascend <- sort(betatilde_N); rank_betatilde_N <- rank(betatilde_N) # (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})'
  delta_beta <- delta_3; changepoints_betatilde <- fc_changepoints(delta_beta, 1, N, betatilde_N_ascend); num_changepoints_betatilde <- length(changepoints_betatilde) # (\widehat{h}_{3,1}, ..., \widehat{h}_{3,\widehat{g}_1})' and \widehat{g}_3 
  ghat_3 <- num_changepoints_betatilde; changepoints_3 <- changepoints_betatilde; rank_3 <- rank_betatilde_N
  ## final estimates (training set)
  thetahat_beta <- fr_thetahat_beta_tau(tau, Y_N_T_train, w_N_T_train, phitilde_N_3_train, N, N, ghat_3, c(1:N), c(1:N), changepoints_3, rank_omegatilde_N, rank_alphatilde_N, rank_3)
  thetahat <- c(omegatilde_N_ascend, alphatilde_N_ascend, thetahat_beta)
  ## CV (validation set)
  w_N_T_test <- fr_weights_NT(Y_N_T_test)
  prediction_error <- fr_loss_panelQR(tau, thetahat, Y_N_T_test, w_N_T_test, N, N, ghat_3, c(1:N), c(1:N), changepoints_3, rank_omegatilde_N, rank_alphatilde_N, rank_3)
  return(prediction_error)
}


