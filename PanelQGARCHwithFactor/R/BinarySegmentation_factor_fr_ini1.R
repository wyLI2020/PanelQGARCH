### selection of thresholds \delta_1, \delta_2, \delta_3
fr_phi_l_N_homo <- function(l, theta_vec, g_vec, changepoints_vec, rank_mat) {
  # transform the parameter vector (\omega_{(1)}, ..., \omega_{(g_1)})', (\alpha_{(1)}, ..., \alpha_{(g_2)})', (\beta_{(1)}, ..., \beta_{(g_3)})', or (\lambda_{l(1)}, ..., \lambda_{l(g_{3+r})})' 
  # to the N \times 1 parameter vector (\omega_{1}, ..., \omega_{N})', (\alpha_{1}, ..., \alpha_{N})', (\beta_{1}, ..., \beta_{N})', or (\lambda_{l1}, ..., \lambda_{lN})' under the homogeneity structure
  g_l <- g_vec[l] 
  rank_l <- rank_mat[,l]
  if(l == 1) {
    changepoints_l <- changepoints_vec[1:g_l]
    phi_l_gl_ascend <- theta_vec[1:g_l] # (\omega_{(1)}, ..., \omega_{(g_1)})', (\alpha_{(1)}, ..., \alpha_{(g_2)})', (\beta_{(1)}, ..., \beta_{(g_3)})', or (\lambda_{l(1)}, ..., \lambda_{l(g_{3+r})})' 
  } else {
    changepoints_l <- changepoints_vec[sum(g_vec[1:(l-1)])+(1:g_l)]
    phi_l_gl_ascend <- theta_vec[sum(g_vec[1:(l-1)])+(1:g_l)] # (\omega_{(1)}, ..., \omega_{(g_1)})', (\alpha_{(1)}, ..., \alpha_{(g_2)})', (\beta_{(1)}, ..., \beta_{(g_3)})', or (\lambda_{l(1)}, ..., \lambda_{l(g_{3+r})})' 
  }
  phi_l_N_ascend <- rep(phi_l_gl_ascend, diff(c(0, changepoints_l))) # (\omega_{1}, ..., \omega_{N})', (\alpha_{1}, ..., \alpha_{N})', (\beta_{1}, ..., \beta_{N})', or (\lambda_{l1}, ..., \lambda_{lN})'
  phi_l_N <- phi_l_N_ascend[rank_l]
  return(phi_l_N)
}
# fr_phi_N_homo <- function(theta_vec, g_vec, changepoints_vec, rank_mat) {
#   # transform the parameter vector \bm\theta = (\omega_{(1)}, ..., \omega_{(g_1)}, \alpha_{(1)}, ..., \alpha_{(g_2)}, \beta_{(1)}, ..., \beta_{(g_3)}, \lambda_{1(1)}, ..., \lambda_{1(g_4)}, ..., \lambda_{r(1)}, ..., \lambda_{r(g_{3+r})})' 
#   # to the N \times (3+r) parameter matrix (\bm{\phi}_1, ..., \bm{\phi}_N)' under the homogeneity structure
#   phi_N_3plusr <- sapply(1:length(g_vec), fr_phi_l_N_homo, theta_vec=theta_vec, g_vec=g_vec, changepoints_vec=changepoints_vec, rank_mat=rank_mat)
#   return(phi_N_3plusr)
# }
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

fr_loss_panelQR <- function(tau, theta_vec, Y_N_T, Factors_r_T, w_N_T, g_vec, changepoints_vec, rank_mat) {
  # the loss function of the final estimator 
  phi_N_3plusr <- fc_phi_N_homo(theta_vec, g_vec, changepoints_vec, rank_mat)
  loss_PanelQR <- fc_loss_panelQR(tau, phi_N_3plusr, Y_N_T, Factors_r_T, w_N_T)
  return(loss_PanelQR)
}

fr2_loss_panelQR <- function(tau, theta_l, theta_1TOlminus1, theta_lplus1TOrplus3, Y_N_T, Factors_r_T, w_N_T, g_vec, changepoints_vec, rank_mat) {
  # the loss function of the final estimator 
  theta_vec <- c(theta_1TOlminus1, theta_l, theta_lplus1TOrplus3)
  phi_N_3plusr <- fc_phi_N_homo(theta_vec, g_vec, changepoints_vec, rank_mat)
  loss_PanelQR <- fc_loss_panelQR(tau, phi_N_3plusr, Y_N_T, Factors_r_T, w_N_T)
  return(loss_PanelQR)
}

fr_theta_l_initial_IN_optim <- function(l, phiini_N_3plusr, g_vec, changepoints_vec) {
  # initial value in optim() for estimating \widehat\bm\omega, \widehat\bm\alpha, \widehat\bm\beta, or \widehat\bm\lambda_l
  phiini_l_N <- phiini_N_3plusr[,l] # (\widetilde\omega_1, ..., \widetilde\omega_N)', (\widetilde\alpha_1, ..., \widetilde\alpha_N)', (\widetilde\beta_1, ..., \widetilde\beta_N)', or (\widetilde\lambda_{lr1}, ..., \widetilde\lambda_{lN})'
  phiini_l_N_ascend <- sort(phiini_l_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})', (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})', (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})', or (\widetilde\lambda_{l(1)}, ..., \widetilde\lambda_{l(N)})'
  g_l <- g_vec[l]
  if(l == 1) {
    changepoints_l <- changepoints_vec[1:g_l]
  } else {
    changepoints_l <- changepoints_vec[sum(g_vec[1:(l-1)])+(1:g_l)]
  }
  theta_l_ini <- phiini_l_N_ascend[changepoints_l]
  return(theta_l_ini)
}
# fr_theta_l_initial_IN_optim <- function(l, phiini_N_3plusr, g_vec, changepoints_vec) {
#   # initial value in optim() for estimating \widehat\bm\omega, \widehat\bm\alpha, \widehat\bm\beta, or \widehat\bm\lambda_l
#   phiini_l_N <- phiini_N_3plusr[,l] # (\widetilde\omega_1, ..., \widetilde\omega_N)', (\widetilde\alpha_1, ..., \widetilde\alpha_N)', (\widetilde\beta_1, ..., \widetilde\beta_N)', or (\widetilde\lambda_{lr1}, ..., \widetilde\lambda_{lN})'
#   phiini_l_N_ascend <- sort(phiini_l_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})', (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})', (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})', or (\widetilde\lambda_{l(1)}, ..., \widetilde\lambda_{l(N)})'
#   g_l <- g_vec[l] 
#   if(l == 1) {
#     changepoints_l <- changepoints_vec[1:g_l]
#   } else {
#     changepoints_l <- changepoints_vec[sum(g_vec[1:(l-1)])+(1:g_l)]
#   }
#   theta_l_ini <- fc_calculate_segment_means(phiini_l_N_ascend, c(0, changepoints_l))
#   return(theta_l_ini)
# }

fr_thetahat_tau <- function(tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, g_vec, changepoints_vec, rank_mat) {
  # the final estimator \widehat\bm\theta = (\widehat\omega_{(1)}, ..., \widehat\omega_{(\widehat{g}_1)}, \widehat\alpha_{(1)}, ..., \widehat\alpha_{(\widehat{g}_2)}, \widehat\beta_{(1)}, ..., \widehat\beta_{(\widehat{g}_3)})', which is obtained using the threshold \delta
  theta_ini <- fc_theta_initial_IN_optim(phitilde_N_3plusr, g_vec, changepoints_vec)
  fit_thetahat <- optim(par=theta_ini, fn=fr_loss_panelQR, method="L-BFGS-B", lower=c(rep(-Inf,g_vec[1]+g_vec[2]),rep(0,g_vec[3]),rep(-Inf,sum(g_vec[-c(1:3)]))), upper=c(rep(Inf,g_vec[1]+g_vec[2]),rep(1-1e-3,g_vec[3]),rep(Inf,sum(g_vec[-c(1:3)]))), tau=tau, Y_N_T=Y_N_T, Factors_r_T=Factors_r_T, w_N_T=w_N_T, g_vec=g_vec, changepoints_vec=changepoints_vec, rank_mat=rank_mat) 
  thetahat <- fit_thetahat$par
  return(thetahat)
}

fr_thetahat_l_tau <- function(l, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, g_vec, changepoints_vec, rank_mat) {
  # if l=1, (\widehat\omega_{(1)}, ..., \widehat\omega_{(\widehat{g}_1)})' when given (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)}, \widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})', \widetilde\lambda_{(11)}, ..., \widetilde\beta_{(1N)})', ..., \widetilde\lambda_{(r1)}, ..., \widetilde\lambda_{(rN)})'
  # if l=2, (\widehat\alpha_{(1)}, ..., \widehat\alpha_{(\widehat{g}_2)})' when given (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)}, \widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})', \widetilde\lambda_{(11)}, ..., \widetilde\beta_{(1N)})', ..., \widetilde\lambda_{(r1)}, ..., \widetilde\lambda_{(rN)})'
  # if l=3, (\widehat\beta_{(1)}, ..., \widehat\beta_{(\widehat{g}_3)})' when given (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)}, \widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})', \widetilde\lambda_{(11)}, ..., \widetilde\beta_{(1N)})', ..., \widetilde\lambda_{(r1)}, ..., \widetilde\lambda_{(rN)})'
  # if l>3, (\widehat\lambda_{l-3,(1)}, ..., \widehat\lambda_{l-3,(\widehat{g}_3)})' when given (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)}, \widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})', \widetilde\lambda_{(11)}, ..., \widetilde\beta_{(1N)})', ..., \widetilde\lambda_{(l-4,1)}, ..., \widetilde\lambda_{(l-4,N)})', ..., \widetilde\lambda_{(l-2,1)}, ..., \widetilde\lambda_{(l-2,N)})', ..., \widetilde\lambda_{(r1)}, ..., \widetilde\lambda_{(rN)})'
  rplus3 <- ncol(phitilde_N_3plusr)
  g_l <- g_vec[l]
  theta_ini <- fc_theta_initial_IN_optim(phitilde_N_3plusr, g_vec, changepoints_vec)
  if(l == 1) {
    theta_ini_l <- theta_ini[1:g_l]
    theta_ini_1TOlminus1 <- vector(mode = "numeric")
    theta_ini_lplus1TOrplus3 <- theta_ini[(g_l+1):length(theta_ini)]
  } else if(l == rplus3) {
    theta_ini_l <- theta_ini[sum(g_vec[1:(l-1)])+(1:g_l)]
    theta_ini_1TOlminus1 <- theta_ini[1:sum(g_vec[1:(l-1)])]
    theta_ini_lplus1TOrplus3 <- vector(mode = "numeric")
  } else {
    theta_ini_l <- theta_ini[sum(g_vec[1:(l-1)])+(1:g_l)]
    theta_ini_1TOlminus1 <- theta_ini[1:sum(g_vec[1:(l-1)])]
    theta_ini_lplus1TOrplus3 <- theta_ini[(sum(g_vec[1:l])+1):length(theta_ini)]
  }
  if(l == 3) {
    lower=rep(0,g_l); upper=rep(1-1e-3,g_l)
  } else {
    lower=rep(-Inf,g_l); upper=rep(Inf,g_l)
  }
  fit_theta_l_hat <- optim(par=theta_ini_l, fn=fr2_loss_panelQR, method="L-BFGS-B", lower=lower, upper=upper, tau=tau, theta_1TOlminus1=theta_ini_1TOlminus1, theta_lplus1TOrplus3=theta_ini_lplus1TOrplus3, Y_N_T=Y_N_T, Factors_r_T=Factors_r_T, w_N_T=w_N_T, g_vec=g_vec, changepoints_vec=changepoints_vec, rank_mat=rank_mat) 
  theta_l_hat <- fit_theta_l_hat$par
  return(theta_l_hat)
}

fr3_BIC_delta_l_tau <- function(l, delta_l, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, rank_mat) {
  # determine the threshold \delta_l using BIC, l = 1, 2, 3, ..., 3+r
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  rplus3 <- ncol(phitilde_N_3plusr)
  ## initial estimates
  phi_l_tilde_N <- phitilde_N_3plusr[,l] # (\widetilde\omega_1, ..., \widetilde\omega_N)', (\widetilde\alpha_1, ..., \widetilde\alpha_N)', (\widetilde\beta_1, ..., \widetilde\beta_N)', or (\widetilde\lambda_{l1}, ..., \widetilde\lambda_{lN})'
  if(l == 1) {
    phi_1TOlminus1_tilde_vec <- vector(mode = "numeric")
    phi_lplus1TOrplus3_tilde_vec <- as.vector(phitilde_N_3plusr[,(l+1):rplus3]) 
  } else if(l == rplus3) {
    phi_1TOlminus1_tilde_vec <- as.vector(phitilde_N_3plusr[,1:(l-1)]) 
    phi_lplus1TOrplus3_tilde_vec <- vector(mode = "numeric")
  } else {
    phi_1TOlminus1_tilde_vec <- as.vector(phitilde_N_3plusr[,1:(l-1)]) 
    phi_lplus1TOrplus3_tilde_vec <- as.vector(phitilde_N_3plusr[,(l+1):rplus3]) 
  }
  ## change points detection
  phi_l_tilde_N_ascend <- sort(phi_l_tilde_N); rank_l <- rank(phi_l_tilde_N) # (\widetilde\omega_{(1)}, ..., \widetilde\omega_{(N)})', (\widetilde\alpha_{(1)}, ..., \widetilde\alpha_{(N)})', (\widetilde\beta_{(1)}, ..., \widetilde\beta_{(N)})', or (\widetilde\lambda_{(l1)}, ..., \widetilde\lambda_{(lN)})'
  changepoints_l <- fc_changepoints(delta_l, 1, N, phi_l_tilde_N_ascend); ghat_l <- length(changepoints_l) # (\widehat{h}_{l,1}, ..., \widehat{h}_{l,\widehat{g}_l})' and \widehat{g}_l 
  ## final estimates
  g_vec <- rep(N, rplus3); g_vec[l] <- ghat_l
  changepoints_1TOlminus1 <- rep(1:N, l-1) 
  changepoints_lplus1TOrplus3 <- rep(1:N, rplus3-l) 
  changepoints_vec <- c(changepoints_1TOlminus1, changepoints_l, changepoints_lplus1TOrplus3)
  thetahat_l <- fr_thetahat_l_tau(l, tau, Y_N_T, Factors_r_T, w_N_T, phitilde_N_3plusr, g_vec, changepoints_vec, rank_mat)
  thetahat <- c(phi_1TOlminus1_tilde_vec, thetahat_l, phi_lplus1TOrplus3_tilde_vec)
  ## BIC
  loss_panelQR <- fr_loss_panelQR(tau, thetahat, Y_N_T, Factors_r_T, w_N_T, g_vec, changepoints_vec, rank_mat)
  # BIC <- 2*N*T*log(loss_panelQR/(N*T)) + length(thetahat)*log(N*T)
  # BIC <- 2*N*log(loss_panelQR/N) + ghat_1*log(N)
  BIC <- loss_panelQR + ghat_l*log(N)
  return(BIC)
}


