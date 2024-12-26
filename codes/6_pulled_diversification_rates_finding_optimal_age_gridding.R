library(castor)
library(ape)
library(dispRity)

# PLEASE CHECK THE FOLLOWING TWO LINES
# if you are using the "Rscript" command in linux/mac terminal
setwd(paste0(getwd(),"/../data/all_trees"))

# if you are using RStudio IDE
#setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../data/all_trees"))

filenames <- list.files()
lineage.names <- substr(filenames, 1, (nchar(filenames)-4))

model_fits <- matrix(nrow=length(filenames),ncol=13)
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

for(i in 1:length(filenames)){
  #i <- 1
  cat('Currently using tree for',lineage.names[i],"\n\n")
  tree <- read.tree(filenames[i])
  # if(Ntip(tree)<8){
  #   cat("Too few tips. Moving to the next tree.\n")
  #   next
  # }
  root_age <- get_tree_span(tree)$max_distance
  
  # time-independent (constant) estimates
  pulled_div_rates_constant <- tryCatch(expr = {suppressWarnings(fit_hbd_pdr_on_grid(tree = tree,
                                                                                     age_grid = NULL,
                                                                                     Ntrials = 5,
                                                                                     Nbootstraps = 5,
                                                                                     Nthreads = 5,
                                                                                     verbose = TRUE))},
                                        error = function(e){return(NA)},
                                        warning = function(w){},
                                        finally = {})
  
  # # time-dependent (varying) estimates
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
  
  pulled.res  <- c()
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
  write.csv(model_fits,"../model_comparisons_parallelized/model_fits_PDR.csv")
}

