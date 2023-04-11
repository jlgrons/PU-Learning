source('SAPUL-functions.R')

run_sims <- function(n_sims = 200, prev = 0.3, n, p, corr=F, col=1){

  all_auc <- c()
  true_beta_all <- c()
  surr_beta_all <- c()
  sass_beta_all <- c()
  sup_100_beta_all <- c()
  sup_200_beta_all <- c()
  sapul_beta_all <- c()
  
  for(i in 1:n_sims){
    
    # Training data.
    my_data <- generate_some_data(p = p, prev = prev, n = n, corr=corr, col=col)
    y <- my_data[, 1]
    s <- my_data[, ncol(my_data)]
    a <- my_data[, 2]
    x <- my_data[, -c(1,2, ncol(my_data))]
    
    # Validation data. 
    my_val_data <- generate_some_data(p = p, prev = prev, n = n, corr=corr, col=col)
    y_v <- my_val_data[, 1]
    s_v <- my_val_data[, ncol(my_data)]
    a_v <- my_val_data[, 2]
    x_v <- my_val_data[, -c(1,2, ncol(my_data))]
    
    # Ideal setting.
    true_beta <- Est.ALASSO.GLM.Approx(cbind(y, s, x))
    true_beta
    
    true_auc <- get_auc(cbind(y_v, cbind(1, s_v, x_v) %*% true_beta))
    true_auc
    
    # Surrogate only setting.
    surr_beta <- Est.ALASSO.GLM(cbind(s, x), regularize = FALSE, fam = 'gaussian')
    surr_beta
    
    surr_auc <- get_auc(cbind(y_v, cbind(1, x_v) %*% surr_beta))
    surr_auc
    
    # Surrogate + labeled setting.
    sass_initial <- glm(y[1:100]~ s[1:100] + x[1:100,] %*% surr_beta[-1], 
                        family = 'binomial')$coef
    sass_beta <- c(sass_initial[c(1,2)], sass_initial[3]*surr_beta[-1])
    
    sass_auc <- get_auc(cbind(y_v, cbind(1, s_v, x_v) %*% sass_beta))
    sass_auc
    
    # Supervised with 100 labeled examples.
    n0 <- 100
    sup_100_beta <- Est.ALASSO.GLM(cbind(y, s, x)[sample(1:length(y), n0), ], regularize = FALSE)
    #Est.ALASSO.GLM.Approx(cbind(y, s, x)[sample(1:length(y), 1000), ])
    sup_100_beta
    
    sup_100_auc <- get_auc(cbind(y_v, cbind(1, s_v, x_v) %*% sup_100_beta))
    sup_100_auc
    
    # Supervised with 200 labeled examples.
    n0 <- 200
    sup_200_beta <- Est.ALASSO.GLM(cbind(y, s, x)[sample(1:length(y), n0), ], regularize = FALSE)
    #Est.ALASSO.GLM.Approx(cbind(y, s, x)[sample(1:length(y), 1000), ])
    sup_200_beta
    
    sup_200_auc <- get_auc(cbind(y_v, cbind(1, s_v, x_v) %*% sup_200_beta))
    sup_200_auc
    
    # Positive + surrogate setting.
    
    ## Step 1: Obtain the linear predictor from the surrogate regression.
    surr_lp <- x %*% surr_beta[-1]
    ## Step 2: Regress the linear predictor + s onto the anchor. 
    initial_param <- c(0.5, 1, 1, 0.5)
    beta_initial <- optim(initial_param, anchor.logitlik.fun, 
                          dat = cbind(a, s, surr_lp))
    beta_initial$par
    
    ## Step 3: Get final estimate. 
    sapul_beta <- c(beta_initial$par[c(2,3)], beta_initial$par[4]*surr_beta[-1])
    sass_initial
    sapul_beta
    sapul_auc <- get_auc(cbind(y_v, cbind(1, s_v, x_v) %*% sapul_beta))
    sapul_auc
    
    # Save everything.
    all_auc <- rbind(all_auc, c(true_auc, surr_auc, sass_auc, sup_100_auc,
                                sup_200_auc, sapul_auc))
    all_auc
    # Save all the betas. 
    true_beta_all <- rbind(true_beta_all, true_beta)
    surr_beta_all <-  rbind(surr_beta_all, true_beta)
    sass_beta_all <-  rbind(sass_beta_all, true_beta)
    sup_beta_100_all <-  rbind(sup_100_beta_all, true_beta)
    sup_beta_200_all <-  rbind(sup_100_beta_all, true_beta)
    sapul_beta_all <-  rbind(sapul_beta_all, true_beta)
    
  }
  
  return(round(colMeans(as.matrix(all_auc)),4))
}

results_matrix <- matrix(0, 8, 6)
ns <- c(100, 200, 50)
ps <- c(50, 100)
for(i in 1:3){
  for(j in 1:2){
    results_matrix[2*i-(j%%2),] = run_sims(n = ns[i], p = ps[j])
  }
}
results_matrix[7,] = run_sims(n = 200, p = 50, corr = T, col = 2)
results_matrix[8,] = run_sims(n = 50, p = 100, corr = T, col = 1)
results <- data.frame(results_matrix)
colnames(results) <- c("true", "surr", "sass", "sup-100", "sup-200","sapul")
results











