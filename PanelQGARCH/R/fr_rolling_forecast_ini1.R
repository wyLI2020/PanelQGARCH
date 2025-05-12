# library(Rcpp)
# library(RcppArmadillo)
# library(parallel)
# source("EstimationProcedure_fr_ini1.R")
# Rcpp::sourceCpp("fc_rolling_forecast.cpp")

fr_ifCover <- function(observation, q) {
  ifCover <- NA
  if(observation < q) {
    ifCover <- 1
  } else {
    ifCover <- 0
  }
  return(ifCover)
}

check_function <- function(x,tau){
  return(x*(tau-(x<=0)))
}

VaRbacktestB <- function(tau, Hit, VaR, Y_true, p){
  # Define a VaR backtest function
  # ============================================================
  # Hit is a sequence of Hit_t=I(y_t<-VaR_t)
  # VaR is a sequence of Value-at-Risk
  # tau is the quantile level of -VaR=Q_tau(y_t|F_{t-1})
  # p is the dimension of lagged hits in DQ tests
  # ============================================================
  n <- length(Hit) # sample size
  ### Unconditional coverage test
  n1 <- sum(Hit); tauhat <- n1/n
  LR_UC <- -2*log((tau/tauhat)^n1*((1-tau)/(1-tauhat))^(n-n1))
  P.UC <- pchisq(LR_UC,df=1,lower.tail=FALSE) # p-value
  ### Independence test
  ContTable <- table(Hit[-n],Hit[-1])
  if(any(dim(ContTable)!=c(2,2))){LR_Ind=NA;P.Ind = NA}else{
    n00 <- ContTable[1,1]; n01 <- ContTable[1,2]; n10 <- ContTable[2,1]; n11 <- ContTable[2,2]
    tau_null <- (n01+n11)/(n00+n01+n10+n11)
    tau0_alt <- n01/(n00+n01); tau1_alt <- n11/(n10+n11)
    LR_Ind <- -2*log((tau_null/tau0_alt)^n01*(tau_null/tau1_alt)^n11*((1-tau_null)/(1-tau0_alt))^n00*((1-tau_null)/(1-tau1_alt))^n10)
    P.Ind <- pchisq(LR_Ind,df=1,lower.tail=FALSE) # p-value
  }
  
  ### Conditional coverage test
  if(any(dim(ContTable)!=c(2,2))){
    LR_CC <- NaN
    P.CC <- NaN
  }else{
    LR_CC <- LR_UC+LR_Ind
    P.CC <- pchisq(LR_CC,df=2,lower.tail=FALSE) # p-value
  }
  ### DQ test: hits
  X <- cbind(1,embed(Hit[-n],p)) # (n-p)*(p+1) matrix
  if (det(t(X)%*%X)<=10^{-5}){P.DQ <- NaN
  } else {
    DQ <- t(Hit[(p+1):n]-tau)%*%X%*%solve(t(X)%*%X)%*%t(X)%*%(Hit[(p+1):n]-tau)/tau/(1-tau)
    P.DQ <- pchisq(DQ,df=p+1,lower.tail=FALSE) # p-value
  }
  
  ECR <- mean(Hit)*100
  PE <- sqrt(n)*abs(mean(Hit)-tau)/sqrt(tau*(1-tau))
  Loss <- mean(check_function(x = (Y_true-(-VaR)),tau = tau))
  ###########################################################################
  # rbind.data.frame(ECR=ECR,PE=PE,Loss=Loss,UC=P.UC,Ind=P.Ind,CC=P.CC,DQ=P.DQ)  # deparse.level=0
  results <- c(ECR, PE, P.CC, P.DQ); names(results) <- c("ECR", "PE", "CC", "DQ")
  return(results)
}

# fr_rolling_forecast <- function(t0, tau, Y_N_T, if_group=c(TRUE, TRUE, TRUE), if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123_candid_lower=c(0,0,0), delta_123_candid_upper=c(2,2,2), testLength=NA, sig_level=0.01, clsNum) {
#   N <- nrow(Y_N_T); T <- ncol(Y_N_T)
#   t1 <- T - t0 # length of the test set
#   SEQ <- 1:t1
#   parFun <- function(l) {
#     Y_N_t0 <- Y_N_T[,(l-1)+(1:t0)]; y_t0plusl <- Y_N_T[,(t0+l)]
#     est_onerep <- tryCatch(fr_onerep_tau(tau, Y_N_t0, if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, testLength), error=function(e){result <- rep(NA, 4*N); return(result)})
#     phitilde_N_3 <- matrix(est_onerep[1:(3*N)], nrow = N, ncol = 3)
#     phihat_N_3 <- matrix(est_onerep[6*N + (1:(3*N))], nrow = N, ncol = 3)
#     qtilde_y_t0plusl <- fc_q_y_t0plus1(Y_N_t0, phitilde_N_3)
#     qhat_y_t0plusl <- fc_q_y_t0plus1(Y_N_t0, phihat_N_3)
#     ifCovers_qtilde_t0plusl <- fc_ifCovers_vec(y_t0plusl, qtilde_y_t0plusl)
#     ifCovers_qhat_t0plusl <- fc_ifCovers_vec(y_t0plusl, qhat_y_t0plusl)
#     VaR_qtilde_t0plusl <- -qtilde_y_t0plusl
#     VaR_qhat_t0plusl <- -qhat_y_t0plusl
#     Hit_qtilde_t0plusl <- ifCovers_qtilde_t0plusl
#     Hit_qhat_t0plusl <- ifCovers_qhat_t0plusl
#     forecast_t0plusl <- c(VaR_qtilde_t0plusl, VaR_qhat_t0plusl, Hit_qtilde_t0plusl, Hit_qhat_t0plusl)
#     return(forecast_t0plusl)
#   }
#   cls <- makeCluster(clsNum, type="FORK")
#   clusterSetRNGStream(cls, 666)
#   forecast_4N_t1 <- parSapply(cls, SEQ, parFun) # Carry out the task parFun parallely and obtain estimates of parameters in Rnum replications
#   stopCluster(cls) # stop workers
#   VaR_qtilde_N_t1 <- forecast_4N_t1[1:N,]
#   VaR_qhat_N_t1 <- forecast_4N_t1[N + (1:N),]
#   Hit_qtilde_N_t1 <- forecast_4N_t1[2*N + (1:N),]
#   Hit_qhat_N_t1 <- forecast_4N_t1[3*N + (1:N),]
#   SEQ_2 <- 1:N
#   parFun_2 <- function(i) {
#     VaR_qtilde_i <- VaR_qtilde_N_t1[i,]
#     VaR_qhat_i <- VaR_qhat_N_t1[i,]
#     Hit_qtilde_i <- Hit_qtilde_N_t1[i,]
#     Hit_qhat_i <- Hit_qhat_N_t1[i,]
#     results_tilde <- VaRbacktestB(tau, Hit_qtilde_i, VaR_qtilde_i, Y_N_T[i,t0+(1:t1)], p=4)
#     results_hat <- VaRbacktestB(tau, Hit_qhat_i, VaR_qhat_i, Y_N_T[i,t0+(1:t1)], p=4)
#     results <- c(results_tilde, results_hat)
#   }
#   cls <- makeCluster(clsNum, type="FORK")
#   clusterSetRNGStream(cls, 666)
#   forecast_summary <- parSapply(cls, SEQ_2, parFun_2) # Carry out the task parFun parallely and obtain estimates of parameters in Rnum replications
#   stopCluster(cls) # stop workers
#   forecast_summary_average <- matrix(NA, nrow = 2, ncol = 4)
#   rownames(forecast_summary_average) <- c("initial", "final")
#   colnames(forecast_summary_average) <- c("ECR", "PE", "CC", "DQ")
#   forecast_summary_average[1,1:2] <- apply(forecast_summary[1:2,], MARGIN = 1, FUN = mean)
#   forecast_summary_average[1,3] <- mean(forecast_summary[3,] > sig_level)
#   forecast_summary_average[1,4] <- mean(forecast_summary[4,] > sig_level)
#   forecast_summary_average[2,1:2] <- apply(forecast_summary[5:6,], MARGIN = 1, FUN = mean)
#   forecast_summary_average[2,3] <- mean(forecast_summary[7,] > sig_level)
#   forecast_summary_average[2,4] <- mean(forecast_summary[8,] > sig_level)
#   # forecast_results <- data.frame("results_individuals"=forecast_summary, "results_average"=forecast_summary_average)
#   return(forecast_summary_average)
# }

fr_rolling_forecast <- function(t0, tau, Y_N_T, if_group=c(TRUE, TRUE, TRUE), if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123_candid_lower=c(0,0,0), delta_123_candid_upper=c(2,2,2), testLength=NA, clsNum) {
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  t1 <- T - t0 # length of the test set
  SEQ <- 1:t1
  parFun <- function(l) {
    Y_N_t0 <- Y_N_T[,(l-1)+(1:t0)]; y_t0plusl <- Y_N_T[,(t0+l)]
    est_onerep <- tryCatch(fr_onerep_tau(tau, Y_N_t0, if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, testLength), error=function(e){result <- rep(NA, 4*N); return(result)})
    phitilde_N_3 <- matrix(est_onerep[1:(3*N)], nrow = N, ncol = 3)
    phihat_N_3 <- matrix(est_onerep[6*N + (1:(3*N))], nrow = N, ncol = 3)
    qtilde_y_t0plusl <- fc_q_y_t0plus1(Y_N_t0, phitilde_N_3)
    qhat_y_t0plusl <- fc_q_y_t0plus1(Y_N_t0, phihat_N_3)
    ifCovers_qtilde_t0plusl <- fc_ifCovers_vec(y_t0plusl, qtilde_y_t0plusl)
    ifCovers_qhat_t0plusl <- fc_ifCovers_vec(y_t0plusl, qhat_y_t0plusl)
    VaR_qtilde_t0plusl <- -qtilde_y_t0plusl
    VaR_qhat_t0plusl <- -qhat_y_t0plusl
    Hit_qtilde_t0plusl <- ifCovers_qtilde_t0plusl
    Hit_qhat_t0plusl <- ifCovers_qhat_t0plusl
    forecast_t0plusl <- c(VaR_qtilde_t0plusl, VaR_qhat_t0plusl, Hit_qtilde_t0plusl, Hit_qhat_t0plusl)
    return(forecast_t0plusl)
  }
  cls <- makeCluster(clsNum, type="FORK")
  clusterSetRNGStream(cls, 666)
  forecast_4N_t1 <- parSapply(cls, SEQ, parFun) # Carry out the task parFun parallely and obtain estimates of parameters in Rnum replications
  stopCluster(cls) # stop workers
  VaR_qtilde_N_t1 <- forecast_4N_t1[1:N,]
  VaR_qhat_N_t1 <- forecast_4N_t1[N + (1:N),]
  Hit_qtilde_N_t1 <- forecast_4N_t1[2*N + (1:N),]
  Hit_qhat_N_t1 <- forecast_4N_t1[3*N + (1:N),]
  SEQ_2 <- 1:N
  parFun_2 <- function(i) {
    VaR_qtilde_i <- VaR_qtilde_N_t1[i,]
    VaR_qhat_i <- VaR_qhat_N_t1[i,]
    Hit_qtilde_i <- Hit_qtilde_N_t1[i,]
    Hit_qhat_i <- Hit_qhat_N_t1[i,]
    results_tilde <- VaRbacktestB(tau, Hit_qtilde_i, VaR_qtilde_i, Y_N_T[i,t0+(1:t1)], p=4)
    results_hat <- VaRbacktestB(tau, Hit_qhat_i, VaR_qhat_i, Y_N_T[i,t0+(1:t1)], p=4)
    results <- c(results_tilde, results_hat)
  }
  cls <- makeCluster(clsNum, type="FORK")
  clusterSetRNGStream(cls, 666)
  forecast_summary <- parSapply(cls, SEQ_2, parFun_2) # Carry out the task parFun parallely and obtain estimates of parameters in Rnum replications
  stopCluster(cls) # stop workers
  return(forecast_summary)
}

fr2_rolling_forecast <- function(t0, tau, Y_N_T, if_group=c(TRUE, TRUE, TRUE), if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123_candid_lower=c(0,0,0), delta_123_candid_upper=c(2,2,2), testLength=NA, clsNum) {
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  t1 <- T - t0 # length of the test set
  SEQ <- 1:t1
  parFun <- function(l) {
    Y_N_t0 <- Y_N_T[,(l-1)+(1:t0)]; y_t0plusl <- Y_N_T[,(t0+l)]
    est_onerep <- tryCatch(fr_onerep_tau(tau, Y_N_t0, if_group, if_mean, criterion, tolerance, num_candid_onesplit, delta_123_candid_lower, delta_123_candid_upper, testLength), error=function(e){result <- rep(NA, 4*N); return(result)})
    phitilde_N_3 <- matrix(est_onerep[1:(3*N)], nrow = N, ncol = 3)
    phihat_N_3 <- matrix(est_onerep[6*N + (1:(3*N))], nrow = N, ncol = 3)
    qtilde_y_t0plusl <- fc_q_y_t0plus1(Y_N_t0, phitilde_N_3)
    qhat_y_t0plusl <- fc_q_y_t0plus1(Y_N_t0, phihat_N_3)
    ifCovers_qtilde_t0plusl <- fc_ifCovers_vec(y_t0plusl, qtilde_y_t0plusl)
    ifCovers_qhat_t0plusl <- fc_ifCovers_vec(y_t0plusl, qhat_y_t0plusl)
    VaR_qtilde_t0plusl <- -qtilde_y_t0plusl
    VaR_qhat_t0plusl <- -qhat_y_t0plusl
    Hit_qtilde_t0plusl <- ifCovers_qtilde_t0plusl
    Hit_qhat_t0plusl <- ifCovers_qhat_t0plusl
    forecast_t0plusl <- c(VaR_qtilde_t0plusl, VaR_qhat_t0plusl, Hit_qtilde_t0plusl, Hit_qhat_t0plusl)
    est_forecast <- c(est_onerep, forecast_t0plusl)
    return(est_forecast)
  }
  cls <- makeCluster(clsNum, type="FORK")
  clusterSetRNGStream(cls, 666)
  est_forecast_13N_t1 <- parSapply(cls, SEQ, parFun) # Carry out the task parFun parallely and obtain estimates of parameters in Rnum replications
  stopCluster(cls) # stop workers
  est_9N_t1 <- est_forecast_13N_t1[(1:(9*N)),]
  forecast_4N_t1 <- est_forecast_13N_t1[-c(1:(9*N)),]
  VaR_qtilde_N_t1 <- forecast_4N_t1[1:N,]
  VaR_qhat_N_t1 <- forecast_4N_t1[N + (1:N),]
  Hit_qtilde_N_t1 <- forecast_4N_t1[2*N + (1:N),]
  Hit_qhat_N_t1 <- forecast_4N_t1[3*N + (1:N),]
  SEQ_2 <- 1:N
  parFun_2 <- function(i) {
    VaR_qtilde_i <- VaR_qtilde_N_t1[i,]
    VaR_qhat_i <- VaR_qhat_N_t1[i,]
    Hit_qtilde_i <- Hit_qtilde_N_t1[i,]
    Hit_qhat_i <- Hit_qhat_N_t1[i,]
    results_tilde <- VaRbacktestB(tau, Hit_qtilde_i, VaR_qtilde_i, Y_N_T[i,t0+(1:t1)], p=4)
    results_hat <- VaRbacktestB(tau, Hit_qhat_i, VaR_qhat_i, Y_N_T[i,t0+(1:t1)], p=4)
    results <- c(results_tilde, results_hat)
  }
  cls <- makeCluster(clsNum, type="FORK")
  clusterSetRNGStream(cls, 666)
  forecast_summary <- parSapply(cls, SEQ_2, parFun_2) # Carry out the task parFun parallely and obtain estimates of parameters in Rnum replications
  stopCluster(cls) # stop workers
  results_list <- list("est_9N_t1"=est_9N_t1, "forecast_4N_t1"=forecast_4N_t1, "forecast_summary"=forecast_summary)
  return(results_list)
}

fr3_rolling_forecast_summary <- function(t0, tau, Y_N_T, forecast_4N_t1, clsNum) {
  N <- nrow(Y_N_T); T <- ncol(Y_N_T)
  t1 <- T - t0 # length of the test set
  VaR_qtilde_N_t1 <- forecast_4N_t1[1:N,]
  VaR_qhat_N_t1 <- forecast_4N_t1[N + (1:N),]
  Hit_qtilde_N_t1 <- forecast_4N_t1[2*N + (1:N),]
  Hit_qhat_N_t1 <- forecast_4N_t1[3*N + (1:N),]
  SEQ_2 <- 1:N
  parFun_2 <- function(i) {
    VaR_qtilde_i <- VaR_qtilde_N_t1[i,]
    VaR_qhat_i <- VaR_qhat_N_t1[i,]
    Hit_qtilde_i <- Hit_qtilde_N_t1[i,]
    Hit_qhat_i <- Hit_qhat_N_t1[i,]
    results_tilde <- VaRbacktestB(tau, Hit_qtilde_i, VaR_qtilde_i, Y_N_T[i,t0+(1:t1)], p=4)
    results_hat <- VaRbacktestB(tau, Hit_qhat_i, VaR_qhat_i, Y_N_T[i,t0+(1:t1)], p=4)
    results <- c(results_tilde, results_hat)
  }
  cls <- makeCluster(clsNum, type="FORK")
  clusterSetRNGStream(cls, 666)
  forecast_summary <- parSapply(cls, SEQ_2, parFun_2) # Carry out the task parFun parallely and obtain estimates of parameters in Rnum replications
  stopCluster(cls) # stop workers
  results_list <- list("forecast_4N_t1"=forecast_4N_t1, "forecast_summary"=forecast_summary)
  return(results_list)
}
