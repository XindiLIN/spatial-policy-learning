library(fields)
library(sp)
library(gstat)
library(GpGp)

#=============================================================
# Function: data_generation_general
#
# Purpose:
#   Generate a dataset with customizable treatment and outcome 
#   distributions for simulation studies.
#
# Arguments:
#   n             : Integer. Number of observations to generate.
#   covariate_gen : Function. A function that generates a data frame
#                   of covariates (size n). User-defined.
#   trt_gen       : Function. A function that takes covariates as input 
#                   and returns treatment assignments (continuous or binary).
#   outcome_gen   : Function. A function that takes covariates and treatment 
#                   as input and returns the outcome values.
#   seed          : Integer or NULL. Optional seed for reproducibility.
#
# Returns:
#   data.frame with columns:
#     - X... : Covariates
#     - trt  : Treatment variable
#     - y    : Outcome variable
#
# Notes:
#   - This is a flexible simulation framework.
#   - You can modify trt_gen and outcome_gen to reflect different 
#     functional forms, noise levels, or heterogeneity.
#=============================================================

data_generation_general <- function(d = 50, 
                                    treatment_mean,      # function to generate treatment vector mean
                                    treatment_fun,       # function to generate treatment vector
                                    outcome_mean,        # function to generate outcome given treatment and covariates mean
                                    outcome_fun,         # function to generate outcome given treatment and covariates
                                    gp_params_U = c(1,0.2,0.5,0),
                                    gp_model_U = "matern_isotropic",
                                    spatial_dependence = TRUE,
                                    spatial_confounding = FALSE,
                                    seed = 1998) {
  
  set.seed(seed)
  # --- 1. Create grid ---
  grid_range <- 1
  grid_x <- seq(0, grid_range, length.out = d)
  grid_y <- seq(0, grid_range, length.out = d)
  grd <- expand.grid(coord_x = grid_x, coord_y = grid_y)
  
  # --- 2. Generate covariates ---
  X <- covariate_fun(grd)
  if(spatial_dependence){
    U <- fast_Gp_sim(gp_params_U, gp_model_U,  locs = grd, 30 )  
  } else {
    U <- rep(0, nrow(grd))
  }
  
  
  
  # --- 3. Generate treatment using user-supplied function ---
  if(spatial_confounding){
    mu_T <- treatment_mean(X = X, U = U)  
  } else{
    mu_T <- treatment_mean(X = X, U = NULL)  
  }
  
  T <- treatment_fun(mu_T = mu_T, locs = grd)
  
  # --- 4. Generate regression mean using user-supplied function ---
  mu_Y <- outcome_mean(X = X, T = T, U = U)
  Y <- outcome_fun(mu_Y = mu_Y, locs = grd)
  
  # --- 7. Return data.frame ---
  data <- list(grd = grd, X = X[,c("X1", "X2", "X3", "X4")], U = U,  mu_T = mu_T, T = T, mu_Y = mu_Y, Y = Y)
  return(data)
}



##### Example Generator Functions


# Generate Covariate
# two of them are continuous
# two of them are normally distributed
covariate_fun <- function(grd, gp_params = c(1,0.2,0.5,0)){
  
  # c(variance, range, smoothness, nugget)
  
  # generate continuous 
  X1 <- fast_Gp_sim(gp_params, "matern_isotropic",  locs = grd, 30 )
  X2 <- fast_Gp_sim(gp_params, "matern_isotropic",  locs = grd, 30 )
  # generate categorical
  categories = 1:4
  
  X3 <- fast_Gp_sim(gp_params, "matern_isotropic",  locs = grd, 30 )
  X4 <- fast_Gp_sim(gp_params, "matern_isotropic",  locs = grd, 30 )
  
  X3 <- cut(X3,
           breaks = quantile(X3, probs = seq(0, 1, length.out = length(categories) + 1)),
           labels = categories, include.lowest = TRUE)
  X4 <- cut(X4,
            breaks = quantile(X4, probs = seq(0, 1, length.out = length(categories) + 1)),
            labels = categories, include.lowest = TRUE)
  
  return(data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4))
}


# Normal treatment

treatment_mean <- function(X, U = NULL,
                           beta_cont = c(0.5, -0.3),
                           beta_cat = seq(-1, 1, length.out = 6), # adjust length if # of dummy vars changes
                           cat_design = NULL) {
  # Continuous covariate effect
  mu_cont <- beta_cont[1] * X$X1 + beta_cont[2] * X$X2
  
  # Categorical covariate effect
  if (is.null(cat_design)) {
    cat_design <- model.matrix(~ X$X3 + X$X4)[, -1]
  }
  mu_cat <- as.numeric(cat_design %*% beta_cat)
  
  # True mean (no noise)
  if(is.null(U)){
    mu_T <- mu_cont + mu_cat  
  } else{
    mu_T <- mu_cont + mu_cat  + 0.5 * U
  }
  
  
  return(as.numeric(mu_T))
}

treatment_fun <- function(mu_T,
                          locs,
                          id_sd = 0.5,
                          treatment_range = c(-4,4)) {
  
  n <- nrow(locs)
  # GP noise
  # T_noise <- fast_Gp_sim(gp_params, gp_model, locs = locs, 30)
  
  # i.i.d. noise
  T_noise <- rnorm(n, mean = 0, sd = id_sd)
  
  # Final treatment
  T <- mu_T + T_noise 
  T <- pmax(T, treatment_range[1])
  T <- pmin(T, treatment_range[2])
  
  # truncate the treatment 
  return(as.numeric(T))
}


##############
# Example
##############
# trt = treatment_fun(x,
#                 beta_cont = c(0.5, -0.5),
#                 beta_cat = NULL,
#                 gp_params = c(1, 0.1, 0.5, 1),
#                 locs = grd,
#                 cat_design = NULL)
# trt_mean = treatment_mean(x, beta_cont = trt$beta_cont, beta_cat = trt$beta_cat)
# fields::image.plot(matrix(trt,50,50) )
# fields::image.plot(matrix(trt_mean,50,50) )
##############



# we want the coefficient to be fixed at the outcome_mean model
outcome_mean <- function(X, T, U,
                         beta_cont = c(0.5, -0.3),  # effect of X1, X2
                         beta_treat = -1,           # effect of treatment
                         beta_cat = seq(-1, 1, length.out = 6),     # effect of categorical dummies
                         cat_design = NULL) {
  
  # Continuous covariate effect
  mu_cont <- beta_cont[1] * X$X1 + beta_cont[2] *X$X2
  
  # Categorical covariate effect
  if (is.null(cat_design)) {
    cat_design <- model.matrix(~X$X3 + X$X4)[, -1]
  }
  
  mu_cat <- as.numeric(cat_design %*% beta_cat)
  
  # True mean outcome
  mu_Y <- beta_treat * T + mu_cont + mu_cat + U
  
  return(mu_Y)
}

## this is for linear outcome

outcome_mean_derivative <- function(X, T, U,
                                    beta_cont = c(0.5, -0.3),
                                    beta_treat = -1,
                                    beta_cat = seq(-1, 1, length.out = 6),
                                    cat_design = NULL) {
  # Derivative wrt T is constant = beta_treat
  rep(beta_treat, length(T))
}


# X is a four-dimensional covariate, first two dimensional are continuous, the second two are categorical
outcome_mean_nonlinear <- function(X, T, U,
                                   beta_cont = c(0.5, -0.3),
                                   beta_treat = -1,           # effect of treatment
                                   beta_cat = seq(-1, 1, length.out = 6),     # effect of categorical dummies
                                   cat_design = NULL){
  
  # Continuous covariate effect
  mu_cont <- beta_cont[1] * X$X1 + beta_cont[2] *X$X2 + X$X2 * X$X1 * (as.numeric(X$X3) - 1)/2
  
  # Categorical covariate effect
  if (is.null(cat_design)) {
    cat_design <- model.matrix(~X$X3 + X$X4)[, -1]
  }
  
  mu_cat <- as.numeric(cat_design %*% beta_cat)
  
  # heterogenous treatment effect
  T_hetero <- T * (0.5 + 0.1 * as.numeric(X$X3) + 0.1 * as.numeric(X$X4)) * (-1)
  
  mu_Y <- T_hetero + mu_cont + mu_cat + U
  
  return(mu_Y)
  
}

outcome_mean_nonlinear_derivative <- function(X, T, U,
                                    beta_cont = c(0.5, -0.3),
                                    beta_treat = -1,
                                    beta_cat = seq(-1, 1, length.out = 6),
                                    cat_design = NULL) {
  # Derivative wrt T is constant = beta_treat
  return((0.5 + 0.1 * as.numeric(X$X3) + 0.1 * as.numeric(X$X4)) * (-1))
}



outcome_fun <- function(mu_Y,
                        locs,
                        id_sd = 1) {
  n <- nrow(locs)
  
  # i.i.d. noise
  Y_noise <- rnorm(n, mean = 0, sd = id_sd)
  
  # Observed outcome
  Y <- mu_Y + Y_noise 
  
  return(as.numeric(Y))
}





