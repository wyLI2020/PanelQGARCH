# Check function
#
# @param x Vector.
# @param tau Double.
#
# @return Vector. \eqn{x\[\tau-I(x<0)\]}
#' @export
#
check_function <- function(x,tau){
  return(x*(tau-(x<=0)))
}

# Subgradient of quantile check function
#
# @param x Double.
# @param tau Double.
#
# @return Double tau-I(x<=0).
#' @export
psi_function <- function(x,tau){
  return(tau-(x<=0))
}



# Compute the bandwidth (i.e. \eqn{h}) for QR
#
# @param n Int. Sample size.
# @param tau Double. Single quantile.
#
# @return Vector. Two bandwidths \eqn{h_{HS}} and \eqn{h_{B}}, i.e. \eqn{(h_{HS},h_{B})}.
#' @export
#'
h_compute <- function(n,tau){
  
  h_B <- n^(-1/5)*((4.5*dnorm(qnorm(tau))^4)/(2*qnorm(tau)^2+1)^2)^(1/5)
  h_HS <- n^(-1/3)*qnorm(1-0.05/2)^(2/3)*((1.5*dnorm(qnorm(tau))^2)/(2*qnorm(tau)^2+1)^2)^(1/3)
  return(c(h_B,h_HS))
}


#' #' Conditional quantile function of \eqn{Y_t}
#' #'
#' #' Conditional quantile function of \eqn{Y_t}, i.e. \eqn{Q_\tau(y_t|\mathcal{F}_{t-1})}.
#' #'
#' #' @param parm Vector. Parameter vector (3-dim), i.e. \eqn{\left(\omega(\tau),\alpha_1(\tau),\beta_1(\tau)\right)^\prime}.
#' #' @param Y Vector. Data.
#' #'
#' #' @return Vector. \eqn{\left(Q_\tau(y_2|\mathcal{F}_{0}),\ldots,Q_\tau(y_{N+1}|\mathcal{F}_{N})\right)^\prime}.
#' #' @export
#' #'
#' q_y <- function(parm,Y){
#'   N <- length(Y)
#'   q <- arch_fft(cst = parm[1],epsilon = Y[1:(N)],lambda = parm[2]*parm[3]^(1:(N)-1)) # length N
#'   return(q)
#' }
#' 
#' #' Transformation function in CQR
#' #'
#' #' @param parm_CQR Vector. Parameter vector (4-dim), i.e. \eqn{(a_0,a_1,b_1,\lambda)^\prime}.
#' #' @param tau Vector. Observations.
#' #'
#' #' @return Vector. Transformed parameter vector (3-dim), i.e. \eqn{(a_0Q_{\lambda}(\lambda)/(1-b_1),a_1Q_{\lambda},b_1)^\prime}.
#' #' @export
#' #'
#' g_tau <- function(parm_CQR,tau){
#'   parm_CQR_g <- array(dim = c(3))
#'   parm_CQR_g[1] <- parm_CQR[1]/(1-parm_CQR[3])*Q_tau_function(tau,parm_CQR[4])
#'   parm_CQR_g[2] <- parm_CQR[2]*Q_tau_function(tau,parm_CQR[4])
#'   parm_CQR_g[3] <- parm_CQR[3]
#'   return(parm_CQR_g)
#' }

# Fast Discrete Fourier Transform (FFT) for ARCH(\eqn{\infty}) models
#
# This function provides a fast algorithm to ARCH(\eqn{\infty}) with methodology in Nielsen and Noël(2020).
#
# @param cst Double. \eqn{\omega(\tau)}.
# @param epsilon Vector. \eqn{(Y_1,Y_2,...,Y_{N-1})}
# @param lambda Vector. \eqn{\alpha_1(\tau)(1,\beta_1(\tau),...,\beta_1^{N-2}(\tau))}, where \eqn{N} is the total number of observations.
#
# @return vector. \eqn{\left(\omega(\tau)+\alpha_1(\tau)\sum_{j=1}^{1}\beta_1^{j-1}(\tau)|Y_{t-j}|,...,\omega(\tau)+\alpha_1(\tau)\sum_{j=1}^{N-1}\beta_1^{j-1}(\tau)|Y_{t-j}|\right)}
#' @export
#
#
# @references Nielsen, M. Ø., and Noël, A. (2021). To infinity and beyond: Efficient computation of ARCH(\eqn{\infty}) models. *Journal of Time Series Analysis*, 42, 338-354.
#
arch_fft <- function(cst, epsilon, lambda){
  iT <- length(epsilon)
  np2 <- nextn(2*iT-1, 2)
  sigma2_arch <- fft(fft(c(lambda, rep(0, np2-iT))) * fft(c(abs(epsilon), rep(0, np2-iT))), inverse = T) / np2;
  sigma2_arch <- cst + sigma2_arch[1:iT]
  return(Re(sigma2_arch))
}

# Estimate GARCH model for initial value
#
# This function is used to provide a reasonable initial value for QR based on fitting a GARCH model with normal innovation term. For detail, see \code{fit1_optim}.
#
# @param Y Vector. Data.
# @param fixTau Vector. Multiple quantile levels.
#
# @return A list includes initial values for QR at multiple quantile levels.
#' @export
#'
Init_par <- function(Y,fixTau){
  r <- sqrt(abs(Y))*((Y>=0)-(Y<0))
  spec=ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                  mean.model = list(armaOrder = c(0, 0),include.mean = FALSE),
                  distribution.model = "norm")
  garch<-ugarchfit(spec = spec, data = r, solver = 'hybrid', fit.control = list(stationarity = 1))
  res <- residuals(garch,standardize=T)
  Q_tau <- sapply(fixTau, function(x) quantile((res)^2*((res>0)-(res<=0)),x))
  return(list("GARCH_par" = garch@fit$coef,"Q_tau"=Q_tau))
}



#' Loss function for self-weighted QR estimation
#'
#' @param parm Vector. Parameter vector (3-dim), i.e. \eqn{(\omega(\tau),\alpha_1(\tau),\beta_1(\tau))^\prime}.
#' @param Y Vector. Data.
#' @param w Vector. Self-weights.
#' @param tau Double. Specific quantile level.
#'
#' @return Double. The value of the loss function for self-weighted QR estimation given a parameter vector.
#' @export
#

Loss_QR <- function(parm,Y,Factors_r_N,w,tau){
  N <- length(Y)
  q <- arch_fft(cst = parm[1],
                epsilon = Y[1:(N-1)],
                lambda = parm[2]*parm[3]^(0:(N-2))) + 
       crossprod(parm[-c(1:3)], Factors_r_N[,2:N])
  l <- check_function(Y[2:N]-q,tau)

  lw <- w[2:N]*l

  return(sum(lw)+w[1]*check_function(Y[1]-parm[1],tau)+crossprod(parm[-c(1:3)], Factors_r_N[,1]))
}



#' Gradient of the loss function for self-weighted QR estimation
#'
#' @param parm Vector. Parameter vector (3-dim), i.e. \eqn{(\omega(\tau),\alpha_1(\tau),\beta_1(\tau))^\prime}.
#' @param Y Vector. Data.
#' @param w Vector. Self-weights.
#' @param tau Double. Specific quantile level.
#'
#' @return Vector. The vector of the gradient of loss function for self-weighted QR estimation given a parameter vector.
#' @export
#

Loss_QR_gr <- function(parm,Y,Factors_r_N,w,tau){

  N <- length(Y)
  beta_c <- parm[3]^((1-1):(N-1-1))
  q <- arch_fft(cst = 0,epsilon = Y[1:(N-1)],lambda = beta_c[1:(N-1)])
  l <- tau-1*(Y[2:N]-parm[1]-parm[2]*q-crossprod(parm[-c(1:3)], Factors_r_N[,2:N]) < 0)
  lw <- -w[2:N]*l
  lw <- c(-w[1]*psi_function(Y[1]-parm[1]-crossprod(parm[-c(1:3)], Factors_r_N[,1]),tau),lw)

  gr_M <- matrix(0,nrow = N,ncol = length(parm))
  gr_M[,1] <- 1

  gr_M[2:N,2] <- q

  beta_2 <- parm[2]*(2:N-1)*parm[3]^(2:N-2)
  gr_M[3:N,3] <- arch_fft(cst = 0,epsilon = Y[1:(N-1)],lambda = beta_2[1:(N-1)])[1:(N-2)]
  
  gr_M[,4:length(parm)] <- t(Factors_r_N)
  
  gr <- c(t(lw)%*%gr_M)

  return(gr)
}

#' An optimization function of self-weighted QR estimation given fixed initial value.
#'
#' This function provide an optimization function of self-weighted QR with a fixed initial value.
#' Either method, derivative optimization using \code{stats::optim()} or derivative-free optimization using \code{dfoptim::hjkb()}, can be used.
#' To avoid error or non-convergence in optimization, we add an innovation to the initial value and re-optimize this problem.
#' If the optimization fails on the fixed initial value,
#' a more complicated optimization function \code{fit1_optim_grid} based on greedy search.
#'
#' @param par Vector. Initial value for optimization.
#' @param Y Vector. Data.
#' @param w Vector. Self-weights.
#' @param tau Double. Specific quantile level.
#' @param lower Vector. Lower bound for parameter. The default value is \code{c(NA,NA,1e-3)}.
#' @param upper Vector. Upper bound for parameter. The default value is \code{c(NA,NA,1-1e-3)}.
#' @param method Character. If \code{method="optim"}, derivative optimization by \code{stats::optim()} is used. If \code{method="dfoptim"}, derivative-free optimization by \code{dfoptim::hjkb()} is used.
#' @param iter_max_1 Int. If the optimization function does not converge or the parameter is at the boundary, then re-optimize. Maximum number of repetitions of this step is \code{iter_max_1}.
#' @param iter_max_2 Int. If the condition in \code{iter_max_1} cannot be satisfied, then relax the boundary condition and estimate again. Maximum number of repetitions of this step is \code{iter_max_2-iter_max_1}.
#' @param seed Double. Random seed is used to generate an innovation to perturb the initial value.
#'
#' @note
#' 1. The selection of initial values.
#'
#' Since the QGARCH model is extended from the classical GARCH model, the initial value can be chosen based on GARCH model. Specifically,
#' \itemize{
#' \item Estimate the parameters of GARCH(1,1) model (2.1) with \eqn{r_t=\sqrt{|y_t|}\left(sgn(y_t)\right)} using Gaussian QMLE and \eqn{Q_\tau(\eta_t)} using empirical quantile of \eqn{\hat{\eta}_t} based on package \code{rugarch}.
#' \item The initial value can be chosen as \eqn{\left(\frac{\hat{a}_0}{1-\hat{b}_1}\hat{Q}_\tau(\varepsilon_t),\hat{a}_1\hat{Q}_\tau(\varepsilon_t),\hat{b}_1\right)}, where \eqn{\hat{\varepsilon}_t=\hat{\eta}_t^2sgn(\hat{\eta}_t)}.
#' }
#'
#' 2. The necessity of random seed.
#'
#' To avoid error or non-convergence in optimization, we add an innovation to the initial value and re-optimize this problem, the random seed is set for reproducible results.
#'
#' @return A list of optimization results returned from the \code{optim()} or \code{hjkb()} function.
#' @export
#'
# fit1_optim <-function(par=NULL,Y,w,tau,lower=c(NA,NA,1e-3),upper=c(NA,NA,1-1e-3),method="optim",iter_max_1=10,iter_max_2=20,seed=1234){
fit1_optim <-function(par=NULL,Y,Factors_r_N,w,tau,method="optim",iter_max_1=10,iter_max_2=20){
  # if(!is.na(seed)){
  #   set.seed(seed)
  # }
  
  r <- nrow(Factors_r_N)
  lower=c(NA,NA,1e-3,rep(NA,r)); upper=c(NA,NA,1-1e-3,rep(NA,r))

  if(is.null(par)){
    init_GARCH <- Init_par(Y = Y,fixTau = tau)
    init_loading <- lm(Y ~ t(Factors_r_N))$coefficients[-1]
    par <- c(init_GARCH$GARCH_par[1]/(1-init_GARCH$GARCH_par[3])*init_GARCH$Q_tau,init_GARCH$GARCH_par[2]*init_GARCH$Q_tau,init_GARCH$GARCH_par[3], init_loading)
  }

  if(method=="dfoptim"){
    fit0 <- try(hjkb(par = par,fn = Loss_QR,Y=Y,Factors_r_N=Factors_r_N,w=w,tau=tau,lower = lower,upper=upper),silent = T)
  }else{
    fit0 <- try(optim(par = par,fn = Loss_QR,gr = Loss_QR_gr,Y=Y,Factors_r_N=Factors_r_N,w=w,tau=tau,method = "L-BFGS-B",lower = lower,upper=upper),silent = T)
  }
  k <- 0
  convergence <- 0
  if("try-error"%in% class(fit0)){
    convergence <- 1
  }else{
    if(fit0$convergence!=0 | any(c(fit0$par[3]<=lower[3],fit0$par[3]>=upper[3]))){
      convergence <- 1
    }
  }
  while (convergence!=0 & k<iter_max_1){
    k <- k+1
    if(method=="dfoptim"){
      fit0 <- try(hjkb(par = par+runif(3+r,-0.1,0.1),fn = Loss_QR,Y=Y,Factors_r_N=Factors_r_N,w=w,tau=tau,lower = lower,upper=upper),silent = T)
    }else{
      fit0 <- try(optim(par = par+runif(3+r,-0.1,0.1),fn = Loss_QR,gr = Loss_QR_gr,Y=Y,Factors_r_N=Factors_r_N,w=w,tau=tau,method = "L-BFGS-B",lower = lower,upper=upper),silent = T)
    }
    convergence <- 0
    if("try-error"%in% class(fit0)){
      convergence <- 1
    }else{
      if(fit0$convergence!=0 | any(c(fit0$par[3]<=lower[3],fit0$par[3]>=upper[3]))){
        convergence <- 1
      }
    }
  }

  while (convergence!=0 & k<iter_max_2){
    k <- k+1
    if(method=="dfoptim"){
      fit0 <- try(hjkb(par = par+runif(3+r,-0.1,0.1),fn = Loss_QR,Y=Y,Factors_r_N=Factors_r_N,w=w,tau=tau,lower = lower,upper=upper),silent = T)
    }else{
      fit0 <- try(optim(par = par+runif(3+r,-0.1,0.1),fn = Loss_QR,Y=Y,Factors_r_N=Factors_r_N,w=w,tau=tau,method = "L-BFGS-B",lower = lower,upper=upper),silent = T)
    }
    convergence <- 0
    if("try-error"%in% class(fit0)){
      convergence <- 1
    }
  }

  if(convergence!=0){fit0 <- fit1_optim_grid(Y = Y,Factors_r_N=Factors_r_N,w = w,tau = tau)}
  return(fit0)
}

fr_phitilde_i_tau <- function(tau, y_i_T, Factors_r_T, w_i_T) {
  fit0 <- fit1_optim(Y=y_i_T,Factors_r_N=Factors_r_T, w=w_i_T, tau=tau)
  phitilde_i <- fit0$par
  return(phitilde_i)
}

#' An optimization function of self-weighted QR given multiple initial values
#'
#' This function extends the optimization function \code{fit1_optim()} given with a single initial value to multiple initial values.
#' The method is based on greedy algorithm. In this framework, we start by listing all the possible parameter values with interval \code{1/D}.
#' Then the loss function is evaluated at all points, and the points with the least loss are selected as the starting points for our optimization.
#'
#'
#' @param Y Vector. Data.
#' @param w Vector. Self-weights.
#' @param tau Double. Specific quantile level.
#' @param lower Vector. Lower bound for parameter. The default value is \code{c(NA,NA,1e-3)}.
#' @param upper Vector. Upper bound for parameter. The default value is \code{c(NA,NA,1-1e-3)}.
#' @param D Int. \code{1/D} is the interval for partition. The default value is 10.
#' @param num_best Int. The number of the selected points. The default value is 10.
#'
#' @return A list of optimization results returned from the \code{optim()} or \code{hjkb()} function.
#' @export
#'
fit1_optim_grid <-function(Y,Factors_r_N,w,tau,D=10,num_best=10){

  r <- nrow(Factors_r_N)
  lower=c(NA,NA,1e-3,rep(NA,r)); upper=c(NA,NA,1-1e-3,rep(NA,r))
  
  beta <- rep(1:(D-1)/D,each = (D-1)^2)
  alpha <- rep(qnorm(rep(1:(D-1)/D,each = (D-1))),times = (D-1))
  omega <- qnorm(rep(1:(D-1)/D,times = (D-1)^2))
  loadings_r <- rep(qnorm(rep(1:(D-1)/D,times = (D-1)^2)), r)
  
  par <- matrix(c(omega,alpha,beta,loadings_r),nrow = (D-1)^3,ncol = 3+r)
  loss_par <- apply(par,1,function(x) Loss_QR(parm=x,Y=Y,Factors_r_N=Factors_r_N,w=w,tau = tau))
  index1 <- order(loss_par)
  if(length(index1)>num_best){index1 <- index1[1:num_best]}

  fun <- function(x){
    k <- 0
    convergence <- 0
    fit0 <- try(optim(par = x,fn = Loss_QR,gr = Loss_QR_gr,Y=Y,Factors_r_N=Factors_r_N,w=w,tau = tau,method = "L-BFGS-B",lower = lower,upper=upper),silent = T)
    if("try-error"%in%class(fit0)==T){
      return(NA)
    }else{
      return(fit0$value)
    }

  }

  l <- apply(par[index1,],1,fun)
  index2 <- which.min(l)

  fit0 <- optim(par = par[index1[index2],],fn = Loss_QR,gr = Loss_QR_gr,Y=Y,Factors_r_N=Factors_r_N,w=w,tau = tau,method = "L-BFGS-B",lower = lower,upper=upper)
  names(fit0$par) <- paste(c("omega","alpha1","beta1",paste0("lambda", 1:r)),"(",tau,")",sep = "")
  return(fit0)
}


