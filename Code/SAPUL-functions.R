# Necessary libraries.

library(mvtnorm)
library(glmnet)
library(glmpath)

# Autoregressive correlation matrix.
autocorr_mat <- function(p = 100, rho = 0.3){
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}

# Expit function.
g.logit <- function(xx){exp(xx)/(exp(xx)+1)}

# Likelihood function.
logitlik.fun <- function(bet.mat, dat, fam = 'binomial'){
  yi = dat[,1]; xi = dat[,-1];
  if(fam == 'binomial'){
    pi.mat = g.logit(cbind(1,xi)%*%bet.mat) ## N x B
    apply(log(pi.mat)*yi + log(1-pi.mat)*(1-yi),2,sum)
  }
}

# Anchor likelihood function.
anchor.logitlik.fun = function(params, dat, fam = 'binomial'){
  yi = dat[,1]; xi = dat[,-1];
  pi.mat = g.logit(cbind(1,xi)%*%params[-1]) 
  c = g.logit(params[1])
  -1*(apply(log(c*pi.mat)*yi + log(1-c*pi.mat)*(1-yi),2,sum))
}

# Computes sums efficiently based on ranks
sum.I <-function(yy, FUN, Yi, Vi=NULL, ties.method = "first") 
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  pos <- rank(c(yy,Yi), ties.method = ties.method)[1:length(yy)]-rank(yy,ties.method = ties.method)
  if (substring(FUN,2,2) == "=") pos <- length(Yi)-pos
  if (!is.null(Vi)) {
    if(substring(FUN,2,2) == "=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

# Repeat vc in dm rows. 
VTM <- function(vc, dm){
  matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}

# Computes the AUC.
get_auc <- function(data){
  dd = data[,1]; xx = data[,2]; 
  n0 = sum(1-dd); n1 = sum(dd)
  x0 = xx[dd == 0]; x1 = xx[dd == 1]
  sum((sum.I(x0, "<=", x1)+sum.I(x0,"<",x1))/2)/(n0*n1)
}

# Function to generate the data.
generate_some_data <- function(n = 100, N = 20000, p = 50, prev = 0.2, 
                               x_cov = 0.4, b0 = NULL, g0 = NULL, corr=F, col=1){
  
  ## n: Number of observed labels (positive)
  ## N: number of unobserved labels
  ## p: number of covariates
  ## prev: phenotype prevalence
  ## x_cov: autocorrelation param to generate covariance matrix
  ## b0: initial value for beta
  ## g0: initial value for gamma
  ## corr: indicator to violate independence assumption.
  ## col: column of design matrix to add to surrogate when corr==T
  
  if(is.null(b0)){
    b0 <- c(-0.6, 0.6, 0.3, -0.3, 0.3, rep(0, p - 5))
    beta <- rbind(1, b0)
  }
  
  if(is.null(g0)){
    g0 <- c(1.5)
    gamma <- rbind(1, g0)
  }
  
  N_tot <- n + N
  
  y <- rbinom(N_tot, 1, prev) # True labels, Bernoulli(prev)
  
  x_mean <- VTM(beta[1, ], N_tot) + y * VTM(beta[2, ], N_tot) # matrix of means for each covariate
  x <- rmvnorm(N_tot, sigma = autocorr_mat(p, x_cov)) + x_mean
  
  s_mean <- VTM(gamma[1, ], N_tot) + y * VTM(gamma[2, ], N_tot) # vector of means for each surrogate
  s <- rnorm(N_tot, mean = s_mean) 
  
  # Violate the assumption that X and S are independent given Y
  if(corr == T){
    s = s + x[,col]
  }
  
  a_ind <- sample(which(y == 1), n) # observed labels
  a <- rep(0, N_tot)
  a[a_ind] <- 1
  
  my_dat <- cbind(y, a, x, s)
  colnames(my_dat) <- c('y', 'a', paste0('x', 1:p), 's')
  return(my_dat)
}


# ------------------------------------ Adaptive LASSO ------------------------------------ #

# function for ALASSO
Est.ALASSO.GLM <- function(data, Wi = NULL, rtn = "EST", nopen.ind = NULL, regularize = yes.regularize,
                          BIC.factor = 0.1, offset = NULL, fam0 = "binomial", 
                          ridge = T){
  
  ## data: n x (p+1) matrix with the first column as the outcome and the remaining p columns as X 
  ## Wi: optional vector of weights
  ## rtn: option to return estimate or entire solution path result
  ## nopen.ind: indices of features to not penalize
  ## regularize: indicator for whether to use ALASSO vs. unpenalized GLM
  ## BIC.factor: factor for modified BIC criterion as in Minnier et al (2011)
  ## offset: vector of offset term
  ## fam0: distributional family
  ## ridge: indicator for whether to use ridge vs. unpenalized GLM in initial fit
  
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn = length(y); pp = ncol(x)
  
  if(is.null(Wi)){Wi=rep(1,nn)}; if(is.null(offset)){offset=rep(0,nn)}
  
  if(regularize){
    
    ## ================================================================================ ##
    ## ridge or unpenalized initial estimator
    ## ================================================================================ ##
    
    lam.ridge = pp/nn; 
    
    if(ridge == T){
      bini = as.vector(coef(glmnet::glmnet(x, y, weights = Wi, alpha = 0, standardize = F, lambda = lam.ridge, 
                                           family = fam0, offset = offset)))}
    else{bini = glm(y~x, weights = Wi, family = fam0, offset = offset)$coeff}
    
    ## ================================================================================ ##
    ## adaptive weights for aLASSO
    ## ================================================================================ ##
    w.b = 1/abs(bini[-1]); x.t = x/VTM(w.b, nrow(x))
    
    ## ================================================================================ ##
    ## glmpath provides solution path for a range of penalty parameters                 ##
    ## ================================================================================ ##
    tmpfit = glmpath::glmpath(x.t, y, nopenalty.subset = nopen.ind, family = fam0, weight = Wi,
                              standardize = F, min.lambda = 0, lambda2 = lam.ridge, offset = offset)
    lam.all = c(seq(min(tmpfit$lambda), max(tmpfit$lambda),length = 500))
    b.all = glmpath::predict.glmpath(tmpfit, s = lam.all, type = "coefficients", mode = "lambda", offset = offset)
    b.all = b.all/VTM(c(1,w.b), nrow(b.all)); m0 = length(lam.all)
    
    ## ================================================================================ ##
    ## calculates degree of freedom for all betas (corresponding to different lam.all)  ##
    ## ================================================================================ ##
    df.all = apply(b.all[,-1,drop=F] != 0, 1, sum); x = as.matrix(x)
    
    ## =============================================================================================================== ##
    ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
    ## =============================================================================================================== ##
    BIC.lam = -2*apply(glmpath::predict.glmpath(tmpfit, newx = x.t,newy = y, s = lam.all, type = "loglik",
                                                mode = "lambda", offset = offset), 2, sum) + 
      min(length(y)^BIC.factor, log(length(y)))*df.all
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  }else{
    ## ========================================================================= ##
    ## if regularize = F, then use standard logistic regression w/o penalization ##
    ## ========================================================================= ##
    bhat = bini = glm(y~x, family = fam0, weight = Wi)$coef; lamhat = 0; lam.all = BIC.lam = b.all = NULL
  }
  
  
  out = c("b" = bhat,"bini" = bini,"lamhat" = lamhat)
  if(rtn=="EST"){return(bhat)}else{
    return(list('bhat' = bhat, 'bini' = bini, 'lamhat' =lamhat,
                "b.all" = b.all,"lam.all" = lam.all,
                "fit" = tmpfit,"BIC.lam" = BIC.lam))}
}


# ------------------------------------ Quadratic Approximation of Likelihood ------------------------------------ #

# function for ALASSO with quadratic approximation to the likelihood
Est.ALASSO.GLM.new <- function(data, Wi = NULL, fit.type, nopen.ind = NULL, BIC.factor = 0.1, 
                              offset = NULL, regularize = T, fam = 'binomial'){
  
  ## data: n x (p+1) matrix with the first column as the outcome and the remaining p columns as X 
  ## Wi: optional vector of weights
  ## fit type: if 'exact' uses standard glmnet, if not uses quadratic approx to likelihood
  ## rtn: option to return estimate or entire solution path result
  ## nopen.ind: indices of features to not penalize
  ## regularize: indicator for whether to use ALASSO vs. unpenalized GLM
  ## BIC.factor: factor for modified BIC criterion as in Minnier et al (2011)
  ## offset: vector of offset term
  ## fam0: distributional family
  ## ridge: indicator for whether to use ridge vs. unpenalized GLM in initial fit
  
  if(fit.type == "exact"){
    Est.ALASSO.GLM(data, Wi = Wi, nopen.ind = nopen.ind, BIC.factor = BIC.factor, offset = offset,
                   fam0 = fam, regularize = regularize, ridge = F)
  }else{
    Est.ALASSO.GLM.Approx(data, Wi = Wi, nopen.ind = nopen.ind, BIC.factor = BIC.factor, offset = offset,
                          fam = fam)
  }
}


Est.ALASSO.GLM.Approx <- function(data, Wi = NULL, rtn = "EST", nopen.ind = NULL, BIC.factor = 0.1,
                                 offset = NULL, fam = 'binomial'){
  
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop = F]; nn = length(y); pp = ncol(x)
  if(is.null(Wi)){Wi = rep(1,nn)}; if(is.null(offset)){offset = rep(0,nn)}
  
  #initial beta
  btilde = glm(y~x, family = fam, weight = Wi)
  
  #Jacobian
  Ahat = solve(summary(btilde)$cov.unscaled); btilde = btilde$coef
  Ahat.half = svd(Ahat); Ahat.half = Ahat.half$u%*%diag(sqrt(Ahat.half$d))%*%t(Ahat.half$u)
  Xtilde = Ahat.half; Ytilde = Ahat.half%*%btilde
  
  #trick to force lasso to estimate adaptive lasso
  w.b = 1/abs(btilde); Xtilde.t = Xtilde/VTM(w.b, nrow(Xtilde))
  
  #lasso
  tmpfit = lars::lars(Xtilde.t, Ytilde, type="lasso", normalize = F, intercept = F)
  lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda), length = 500))
  b.all = lars::predict.lars(tmpfit, s = lam.all, type = "coefficients", mode = "lambda")$coef
  b.all = b.all/VTM(w.b,nrow(b.all)); m0 = length(lam.all)
  df.all = apply(b.all[,-1,drop=F] != 0, 1, sum) + 1;
  
  ## =============================================================================================================== ##
  ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
  ## =============================================================================================================== ##
  BIC.lam = -2*logitlik.fun(t(b.all), dat = data, fam = fam) + min(nn^BIC.factor, log(nn))*df.all
  m.opt = (1:m0)[BIC.lam == min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  bhat
}

