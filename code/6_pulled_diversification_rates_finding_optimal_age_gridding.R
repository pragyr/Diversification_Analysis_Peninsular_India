library(castor)
library(ape)
library(dispRity)

setwd("data/all_trees")


filenames <- list.files() 
lineage.names <- substr(filenames, 1, (nchar(filenames)-4))

model_fits <- matrix(nrow=length(filenames),ncol=13) # blank matrix to store the estimates, AIC and BIC of different models
# setting the column names of the matrix
colnames(model_fits) <- c("Lineage name",
                          "rate_constant_pd_rate",
                          "rate_constant_AIC",
                          "rate_constant_BIC",
                          "best_grid_size",
                          "rate_varying_AIC",
                          "rate_varying_BIC",
                          "rate_constant_ps_rate",
                          "rate_constant_AIC_ps",
                          "rate_constant_BIC_ps",
                          "best_grid_size_ps",
                          "rate_varying_AIC_ps",
                          "rate_varying_BIC_ps")

# looping over all the files
for(i in 1:length(filenames)){
  cat('Currently using tree for',lineage.names[i],"\n\n")
  tree <- read.tree(filenames[i]) # reading a tree from a file
  root_age <- get_tree_span(tree)$max_distance # extracting the crown age
  
  # Pulled diversification rates
  # fitting the time-independent (constant) Pulled-Diversification model and storing the results in the following variable
  pulled_div_rates_constant <- tryCatch(expr = {suppressWarnings(fit_hbd_pdr_on_grid(tree = tree,
                                                                                     age_grid = NULL,
                                                                                     Ntrials = 5,
                                                                                     Nbootstraps = 5,
                                                                                     Nthreads = 5,
                                                                                     verbose = TRUE))},
                                        error = function(e){return(NA)},
                                        warning = function(w){},
                                        finally = {}) # using try catch to avoid breaking the loop when errors get raised for an input tree
  
  # fitting the time-dependent (varying) Pulled-Diversification model and storing the results in the following variable
  pulled_div_rates_time_varying <- tryCatch(expr = {suppressWarnings(fit_hbd_pdr_on_best_grid_size(tree,
                                                                                                   grid_sizes = 1:5,
                                                                                                   criterion             = "BIC",
                                                                                                   splines_degree        = 2,
                                                                                                   max_model_runtime = (30*60),
                                                                                                   Ntrials               = 5,
                                                                                                   Nbootstraps           = 5,
                                                                                                   Nthreads              = 5,
                                                                                                   fit_control           = list(),
                                                                                                   verbose               = TRUE,
                                                                                                   verbose_prefix        = ""))},
                                            error = function(e){return(NA)},
                                            warning = function(w){},
                                            finally = {})
  
  pulled.res  <- c() # make a blank vector to store the relevant info from both the result variable defined above
  # if-else construct to check which analysis could not be run and accordingly fill the vector mentioned in the previous line
  # with "-" (when analysis could not run) or the respective values.
  if(is.na(pulled_div_rates_constant)[[1]] | is.na(pulled_div_rates_time_varying)[[1]]){
    if(!is.na(pulled_div_rates_constant)[[1]]){
      cat("ERROR: Time-varying model could not be run at all")
      pulled.res <- c(lineage.names[i],
                      pulled_div_rates_constant$fitted_PDR,
                      pulled_div_rates_constant$AIC,
                      pulled_div_rates_constant$BIC,
                      "-",
                      "-",
                      "-")
    }else{
      cat("ERROR: Could not be run at all")
      pulled.res <- c(lineage.names[i],
                      "-",
                      "-",
                      "-",
                      "-",
                      "-",
                      "-")
    }
    
  }else if(!pulled_div_rates_time_varying$success){
    cat(sprintf("ERROR: Fitting failed: %s\n",pulled_div_rates_time_varying$error))
    pulled.res <- c(lineage.names[i],
                    pulled_div_rates_constant$fitted_PDR,
                    pulled_div_rates_constant$AIC,
                    pulled_div_rates_constant$BIC,
                    "-",
                    "-",
                    "-")
  }else{
    best_fit = pulled_div_rates_time_varying$best_fit
    cat(sprintf("Fitting succeeded:\nBest grid size=%d\n",length(best_fit$age_grid)))
    pulled.res <- c(lineage.names[i],
                    pulled_div_rates_constant$fitted_PDR,
                    pulled_div_rates_constant$AIC,
                    pulled_div_rates_constant$BIC,
                    length(best_fit$age_grid),
                    best_fit$AIC,
                    best_fit$BIC)
  }
  
  ### Pulled speciation rate
  # time-independent (constant) estimates
  pulled_spe_rates_constant <- tryCatch(expr = {suppressWarnings(fit_hbd_psr_on_grid(tree = tree,
                                                                                     age_grid = NULL,
                                                                                     Ntrials = 5,
                                                                                     Nbootstraps = 5,
                                                                                     Nthreads = 5,
                                                                                     verbose = TRUE))},
                                        error = function(e){return(NA)},
                                        warning = function(w){},
                                        finally = {})
  
  # # time-dependent (varying) estimates
  pulled_spe_rates_time_varying <- tryCatch(expr = {suppressWarnings(fit_hbd_psr_on_best_grid_size(tree,
                                                                                                   grid_sizes = 1:5,
                                                                                                   criterion             = "BIC",
                                                                                                   splines_degree        = 2,
                                                                                                   max_model_runtime = (30*60),
                                                                                                   Ntrials               = 5,
                                                                                                   Nbootstraps           = 5,
                                                                                                   Nthreads              = 5,
                                                                                                   fit_control           = list(),
                                                                                                   verbose               = TRUE,
                                                                                                   verbose_prefix        = ""))},
                                            error = function(e){return(NA)},
                                            warning = function(w){},
                                            finally = {})
  
  
  if(is.na(pulled_spe_rates_constant)[[1]] | is.na(pulled_spe_rates_time_varying)[[1]]){
    if(!is.na(pulled_spe_rates_constant)[[1]]){
      cat("ERROR: Time-varying model could not be run at all")
      pulled.res <- c(pulled.res,
                      pulled_spe_rates_constant$fitted_PSR,
                      pulled_spe_rates_constant$AIC,
                      pulled_spe_rates_constant$BIC,
                      "-",
                      "-",
                      "-")
    }else{
      cat("ERROR: Could not be run at all")
      pulled.res <- c(pulled.res,
                      "-",
                      "-",
                      "-",
                      "-",
                      "-",
                      "-")
    }
    
  }else if(!pulled_spe_rates_time_varying$success){
    cat(sprintf("ERROR: Fitting failed: %s\n",pulled_spe_rates_time_varying$error))
    pulled.res <- c(pulled.res,
                    pulled_spe_rates_constant$fitted_PSR,
                    pulled_spe_rates_constant$AIC,
                    pulled_spe_rates_constant$BIC,
                    "-",
                    "-",
                    "-")
  }else{
    best_fit = pulled_spe_rates_time_varying$best_fit
    cat(sprintf("Fitting succeeded:\nBest grid size=%d\n",length(best_fit$age_grid)))
    pulled.res <- c(pulled.res,
                    pulled_spe_rates_constant$fitted_PSR,
                    pulled_spe_rates_constant$AIC,
                    pulled_spe_rates_constant$BIC,
                    length(best_fit$age_grid),
                    best_fit$AIC,
                    best_fit$BIC)
  }
  
  
  ### pulled speciation rate
  # time-independent (constant) estimates
  pulled_spe_rates_constant <- tryCatch(expr = {suppressWarnings(fit_hbd_psr_on_grid(tree = tree,
                                                                                     age_grid = NULL,
                                                                                     Ntrials = 5,
                                                                                     Nbootstraps = 5,
                                                                                     Nthreads = 5,
                                                                                     verbose = TRUE))},
                                        error = function(e){return(NA)},
                                        warning = function(w){},
                                        finally = {})
  
  # # time-dependent (varying) estimates
  pulled_spe_rates_time_varying <- tryCatch(expr = {suppressWarnings(fit_hbd_psr_on_best_grid_size(tree,
                                                                                                   grid_sizes = 1:5,
                                                                                                   criterion             = "BIC",
                                                                                                   splines_degree        = 2,
                                                                                                   max_model_runtime = (30*60),
                                                                                                   Ntrials               = 5,
                                                                                                   Nbootstraps           = 5,
                                                                                                   Nthreads              = 5,
                                                                                                   fit_control           = list(),
                                                                                                   verbose               = TRUE,
                                                                                                   verbose_prefix        = ""))},
                                            error = function(e){return(NA)},
                                            warning = function(w){},
                                            finally = {})
  
  
  if(is.na(pulled_spe_rates_constant)[[1]] | is.na(pulled_spe_rates_time_varying)[[1]]){
    if(!is.na(pulled_spe_rates_constant)[[1]]){
      cat("ERROR: Time-varying model could not be run at all")
      pulled.res <- c(pulled.res,
                      pulled_spe_rates_constant$fitted_PSR,
                      pulled_spe_rates_constant$AIC,
                      pulled_spe_rates_constant$BIC,
                      "-",
                      "-",
                      "-")
    }else{
      cat("ERROR: Could not be run at all")
      pulled.res <- c(pulled.res,
                      "-",
                      "-",
                      "-",
                      "-",
                      "-",
                      "-")
    }
    
  }else if(!pulled_spe_rates_time_varying$success){
    cat(sprintf("ERROR: Fitting failed: %s\n",pulled_spe_rates_time_varying$error))
    pulled.res <- c(pulled.res,
                    pulled_spe_rates_constant$fitted_PSR,
                    pulled_spe_rates_constant$AIC,
                    pulled_spe_rates_constant$BIC,
                    "-",
                    "-",
                    "-")
  }else{
    best_fit = pulled_spe_rates_time_varying$best_fit
    cat(sprintf("Fitting succeeded:\nBest grid size=%d\n",length(best_fit$age_grid)))
    pulled.res <- c(pulled.res,
                    pulled_spe_rates_constant$fitted_PSR,
                    pulled_spe_rates_constant$AIC,
                    pulled_spe_rates_constant$BIC,
                    length(best_fit$age_grid),
                    best_fit$AIC,
                    best_fit$BIC)
  }
  
  model_fits[i,] <- pulled.res
  
  # there has to be a folder by the name of "model_comparisons_parallelized" to execute the next line
  # that folder gets created during the execution of Code 3
  write.csv(model_fits,"../model_comparisons_parallelized/model_fits_PDR.csv") # keep writing the matrix into a CSV file after every iteration
                                                                               # so that the user can monitor the results while the loop is running
}

# set the working directory back to the root of the R project
setwd("../../")