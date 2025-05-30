library(iterators)
library(parallel)
library(doParallel)
library(foreach)
library(cluster)

library(RPANDA)
library(ape)
library(phytools)
library(DDD)
library(dispRity)

set.seed(20)
setwd("data/all_trees")

#### Necessary : Make a folder "all_trees" and store all the tree files in it that would be analysed,
#                Make a folder "paleoclimate_data" and store all the smoothened paleoclimate data in it

#### Input : All trees in a folder : "all_trees", 
#            All smoothened paleoclimate data : "paleoclimate_data"

############################################## RPANDA  MODELS ####################################################


filenames <- list.files()

lineage.names <- substr(filenames,1,nchar(filenames)-4)

#before going further make a directory (outside your "all_trees" directory) for storing tables of model comparisons
# but check if the folder already exists
if(!("model_comparisons_parallelized" %in% list.files("../"))) system2("mkdir","../model_comparisons_parallelized/")
lineage.name <- lineage.names[1]
main_func <- function(lineage.name){
  # had to reload the packages as within the cluster made for parallelisation, everything outside of it becomes obsolete within it
  library(RPANDA)
  library(ape)
  library(phytools)
  library(DDD)
  library(dispRity)
  
  Temp_dat <- read.csv("../paleoclimate_data/Smooth_Inftemp.csv")[,c(2,3)]
  C_dat_potwar <- read.csv("../paleoclimate_data/Smooth_C_dat_Potwar.csv")[,c(2,3)]
  Him.oro_dat <- read.csv("../paleoclimate_data/Smooth_Him_oro_dat.csv")[,c(2,3)]
  
  # Estimate positions of potential peaks and dips in diversification rates (relevant to the function I designed to model symmetrical episodic rate shifts)
  estimate_peak_or_dip_position <- function(Tr,tot_time,res.const.bd,fix.width,peak.in.which.rate){
    res.epi <- c()
    #fix.width = TRUE;peak.in.which.rate = "birth"
    timeseries <- seq(1,tot_time,0.5)
    if(max(timeseries) < tot_time) timeseries <- c(timeseries,tot_time)
    
    if(fix.width){
      cols.mat <- 2
      y <- vector("numeric",2)
    }else{
      cols.mat <- 3
      y <- vector("numeric",3)
    }
    rate_matrix <- matrix(nrow=length(timeseries),ncol=cols.mat)
    
    for(p in timeseries){
      #p=8
      if(fix.width){
        # peak function
        f.epi <- function(t,r){r[1] + (r[2]/p) * (t/((t-p)^2+1))}
      }else{
        y[3] <- 1
        
        # peak function
        f.epi <- function(t,r){r[1] + (r[2]*r[3]/p) * (t/((t-p)^2+r[3]))}
      }
      
      if(peak.in.which.rate == "birth"){
        cst.lamb = FALSE; cst.mu = TRUE; fix.mu = TRUE;
        y[1] <- abs(res.const.bd$lamb_par)
        y[2] <- runif(1,0,(y[1]/10))
        
        lamb_par <- y
        mu_par <- abs(res.const.bd$mu_par)
        
        f.lambda <- f.epi
        f.mu <- function(t,r){abs(r[1])}
      }else if(peak.in.which.rate == "death"){
        cst.lamb = TRUE;cst.mu = FALSE;fix.mu = FALSE
        y[1] <- abs(res.const.bd$mu_par)
        y[2] <- runif(1,0,(y[1]/10))
        
        lamb_par <- abs(res.const.bd$lamb_par)
        mu_par <- y
        
        f.lambda <-  function(t,r){abs(r[1])}
        f.mu <- f.epi
      }else{
        cat("Please check if you wrote either \"birth\" or \"death\".\nPossibly there was a typo.")
        
        return(as.data.frame(t(as.matrix(c("p"=NA,"y"=c(NA,NA,NA))))))
      }
      
      tryCatch(
        expr =  {
          res.t <- fit_bd(Tr,tot_time,f.lambda,f.mu,lamb_par,mu_par,f=1,dt=1e-3,
                          cst.lamb = cst.lamb,cst.mu = cst.mu,fix.mu = fix.mu, expo.lamb = FALSE,expo.mu = FALSE)
          res.epi <- c(res.epi,res.t$aicc)
          
          if(peak.in.which.rate == "birth") rate_matrix[which(timeseries == p),] <- res.t$lamb_par else rate_matrix[which(timeseries == p),] <- res.t$mu_par
          
        },
        error=function(e){},
        warning=function(w){},
        finally = {}
      )
      
    }
    
    # Further processing to get positions for peaks and dips in rates
    # removing aicc values that are negative
    if(length(which(res.epi<0))>0){
      aicc.epi <- res.epi[-which(res.epi<0)]
      timeseries <- timeseries[-which(res.epi<0)]
      rate_matrix <- rate_matrix[-which(res.epi<0),]
    }else aicc.epi <- res.epi
    
    if(!is.null(aicc.epi)){
      p <- timeseries[which(aicc.epi < (res.const.bd$aicc - 2))]
      if(length(p) == 0){
        p <- NA
        cat("\nLocal peaks or dips aren't possibly there in",peak.in.which.rate,"rate\n")
        y <- c(NA,NA,NA)
      }else{
        p <- timeseries[which(aicc.epi == min(aicc.epi))]
        if(length(p) > 1) p <- mean(p)
        cat("\nPotential position(s) for local peaks or dips in",peak.in.which.rate,"rate:",p," Mya\n
        And paramter estimates are:",rate_matrix[which(timeseries == p),])
        y <- rate_matrix[which(timeseries == p),]
      }
      
    }else {
      p <- NA
      y <- c(NA,NA,NA)
      cat("\nLocal peaks or dips aren't possibly there in",peak.in.which.rate," rate\n")
    }
    return(as.data.frame(t(as.matrix(c("p"=p,"y"=y)))))
  }
  
  # Function to fit and find the best parameter estimates for an environmental-variable-dependent diversification models 
  run_fit_env <- function(data,myTree,lineage.name,env.name,vector.b.funcs,vector.d.funcs,res.const.bd){
    times <- branching.times(myTree)
    tot_time<-max(times)
    sorted.times <- rev(as.numeric(sort(times)))
    
    
    # Let's calculate the initial speciation rate (dummy) - rate at crown
    ini.lamb <- (1/(sorted.times[1]-sorted.times[2]))/3
    # # Let's calculate the final speciation rate (dummy) - rate at tips
    fin.lamb <- (1/(sorted.times[length(sorted.times)-1] - sorted.times[length(sorted.times)]))/Ntip(myTree)
    
    
    ########################### Environmental-variable-dependent models ##############################
    ##################################################################################################
    
    
    cat("\nSetting initial parameters for Environmental_variable(",env.name,")-dependent model.",
        "\nCalculating likelihood for Environmental_variable(",env.name,")-dependent model......")
    
    
    LH.env.diff.funcs <- aicc.env.diff.funcs <- lamb.env <- alpha.env <- mu.env <- beta.env <- complexity <-rep(NA,length(vector.b.funcs)*length(vector.d.funcs))
    res.count <- 1
    models.labels <- c()
    
    # running over all possible environmental-variable-dependent models for "Tr" phylogenetic tree
    for(m in 1:length(vector.b.funcs)){
      #m=2
      
      # setting initial parameters for lambda and algorithm
      if(substr(labels(vector.b.funcs)[m],3,5) == "cst"){
        lamb_par_init.env <- res.const.bd$lamb_par
        cst.lamb <- TRUE
        expo.lamb <- FALSE
      }else{
        lamb_par_init.env <- c(fin.lamb, runif(1,0,fin.lamb/2))
        cst.lamb <- FALSE
        if(substr(labels(vector.b.funcs)[m],3,5) == "exp"){
          expo.lamb <- TRUE
        }else{
          expo.lamb <- FALSE
        }
        
      }
      
      # Environmental_variable-dependent function for lambda in this iteration
      lamb.f <- vector.b.funcs[[m]]
      
      for(n in 1:length(vector.d.funcs)){
        #n=1
        if(m==1){
          if(n %in% 1:2) next
        }
        
        # setting initial parameters for mu and algorithm
        if(substr(labels(vector.d.funcs)[n],3,3) == "0"){
          mu_par_init.env <- c()
          fix.mu <- TRUE
          cst.mu <- TRUE
          expo.mu <- FALSE
        }else if(substr(labels(vector.d.funcs)[n],3,5) == "cst"){
          mu_par_init.env <- res.const.bd$mu_par
          cst.mu <- TRUE
          fix.mu <- FALSE
          expo.mu <- FALSE
        }else{
          mu_par_init.env <- c(res.const.bd$mu_par, runif(1,0,0.01))
          cst.mu <- FALSE
          fix.mu <- FALSE
          if(substr(labels(vector.d.funcs)[n],3,5) == "exp"){
            expo.mu <- TRUE
          }else{
            expo.mu <- FALSE
          }
        }
        
        
        # Environmental_variable-dependent function for mu in this iteration
        mu.f <- vector.d.funcs[[n]]
        
        error_counter <- 0      
        tryCatch(
          expr = {
            res.env <- fit_env(myTree,data,max(branching.times(myTree)),lamb.f,mu.f,lamb_par_init.env,mu_par_init.env,f=1,dt=1e-3,
                               cst.lamb = cst.lamb,cst.mu = cst.mu,fix.mu = fix.mu, expo.lamb = expo.lamb,expo.mu = expo.mu)
            
            complexity[res.count] <- length(res.env$lamb_par) + length(res.env$mu_par)
            lamb.env[res.count] <- res.env$lamb_par[1]
            if(substr(labels(vector.b.funcs)[m],3,5) != "cst"){
              alpha.env[res.count] <- res.env$lamb_par[2]
            }
            if(substr(labels(vector.d.funcs)[n],3,3) != "0"){
              mu.env[res.count] <- res.env$mu_par[1]
              
              if(substr(labels(vector.d.funcs)[n],3,5) != "cst"){
                beta.env[res.count] <- res.env$mu_par[2]
              }
            }
            
            LH.env.diff.funcs[res.count] <- res.env$LH
            aicc.env.diff.funcs[res.count] <- res.env$aicc
            
          },
          error = function(e){
            #cat("\nError   ",error_counter,"\n")
          },
          warning = function(w){},
          finally = {}
        )
        
        models.labels <- c(models.labels, paste0(labels(vector.b.funcs)[m],"_",labels(vector.d.funcs)[n]))
        cat("\nLikelihood for ",models.labels[res.count],", of complexity - ",complexity[res.count]," is:",LH.env.diff.funcs[res.count]," and AICc is:",aicc.env.diff.funcs[res.count],"\n")
        res.count = res.count + 1
      }
    }
    
    # choosing best model using information about AICc scores and model complexity
    if(length(which(!is.na(aicc.env.diff.funcs))) != 0){
      min.aicc.env <- min(aicc.env.diff.funcs[!is.na(aicc.env.diff.funcs)]) # get the minimum AICc score
      cat("\nLeast AICc score:",min.aicc.env,"\n")
      
      # check if there are other models that aren't significantly worse than the one with lowest AICc
      ind.equi.models <- which(aicc.env.diff.funcs <= (min.aicc.env + 2)) # indices of equivalent models
      ind.best.model <- ind.equi.models
      
      eq.models.complexity <- complexity[ind.equi.models]
      eq.models.labels <- models.labels[ind.equi.models]
      eq.models.aicc <- aicc.env.diff.funcs[ind.equi.models]
      
      if(length(eq.models.labels)>1){
        cat("\nThere are multiple models with similar AICc scores.\n
          Complexity of best equivalent models -",eq.models.labels,":",eq.models.complexity,"respectively\n\n")
        
        ind.least.complexity <- which(models.labels %in% eq.models.labels[which(eq.models.complexity == min(eq.models.complexity))])
        ind.best.model <- ind.least.complexity
        
        # if complexities are also same then this section gets implemented
        if(length(ind.least.complexity)>1){
          #check if AICc values also exactly match
          ind.least.aicc <- which(models.labels %in% eq.models.labels[which(eq.models.aicc == min.aicc.env)])
          ind.best.model <- ind.least.aicc
          
          if(length(ind.least.aicc)>1){ # if AICc values of multiple models match exactly, then randomly choose one of the models, as all of them are equally likely
            ind.least.aicc <- ind.least.complexity[as.integer(runif(n = 1,min = 1,max = length(ind.least.complexity)))]
            ind.best.model <- ind.least.aicc
          }
        }
      }
      
      min.aicc.env <- aicc.env.diff.funcs[ind.best.model]
      max.LH.env <- LH.env.diff.funcs[ind.best.model]
      max.lamb.env <- lamb.env[ind.best.model]
      max.alpha.env <- alpha.env[ind.best.model]
      max.mu.env <- mu.env[ind.best.model]
      max.beta.env <- beta.env[ind.best.model]
      best.model.env <- models.labels[ind.best.model]
      
      cat("\nBest-fitting model:",best.model.env,"\n")
      res.env <- as.data.frame(t(c(best.model.env,max.lamb.env,max.alpha.env,max.mu.env,max.beta.env,min.aicc.env,max.LH.env)))
      colnames(res.env) <- c("Best Env model","lambda","alpha","mu","beta","aicc","LH")
    }else{
      res.env <- as.data.frame(t(c("Best Env model"=NA,"lambda"=NA,"alpha"=NA,"mu"=NA,"beta"=NA,"aicc"=NA,"LH"=NA)))
    }
    return(res.env)
  }
  
  
  # Function to fit and find the best parameter estimates for an time-varying and diversity-dependent diversification models  
  run_model_fit_wo_env <- function(Tr){
    #Tr <- myTree
    times <- branching.times(Tr)
    tot_time <- max(times)
    
    sorted.times <- rev(as.numeric(sort(times)))
    
    # Let's calculate the initial speciation rate (dummy) - rate at crown
    ini.lamb <- (1/(sorted.times[1]-sorted.times[2]))/3
    # # Let's calculate the final speciation rate (dummy) - rate at tips
    fin.lamb <- (1/(sorted.times[length(sorted.times)-1] - sorted.times[length(sorted.times)]))/Ntip(Tr)
    
    ###################### Constant rate only birth model
    f.lamb.const.b <- function(t,y){abs(y[1])}
    f.mu.const.b <- function(t,y){0}
    
    lamb_par_init.const.b <- dlnorm(runif(1,0,1),meanlog = 0, sdlog = 1, log = FALSE)
    mu_par_init.const.b <- c() 
    
    cat("\nInitial parameters set for CB model.\nCalculating likelihood for CB model......\n")
    
    res.const.b <- fit_bd(Tr,max(branching.times(Tr)),f.lamb.const.b,f.mu.const.b,lamb_par_init.const.b,mu_par_init.const.b,f=1,cst.lamb = TRUE,fix.mu = TRUE)
    cat("\nLikelihood of CB model calculated\n")
    cat("\nEstimated Lambda is: ",res.const.b$lamb_par,"\nMean of calculated initial rate and final rate is: ",mean(c(ini.lamb,fin.lamb)))
    
    ###################### Constant rate Birth and Death model
    f.lamb.const.bd <- function(t,y){abs(y[1])}
    f.mu.const.bd <- function(t,y){abs(y[1])}
    
    lamb_par_init.const.bd <- res.const.b$lamb_par
    mu_par_init.const.bd <-  runif(1,0,1)
    while(mu_par_init.const.bd >= lamb_par_init.const.bd){
      mu_par_init.const.bd <-  runif(1,0,1)
    }
    
    cat("\nInitial parameters set for CBD model.\nCalculating likelihood for CBD model......\n")
    
    res.const.bd <- fit_bd(Tr,max(branching.times(Tr)),f.lamb.const.bd,f.mu.const.bd,lamb_par_init.const.bd,mu_par_init.const.bd,f=1,cst.lamb = TRUE,cst.mu = TRUE)
    cat("\nLikelihood of CBD model calculated\n")
    
    cat("\nEstimated Lambda is: ",res.const.bd$lamb_par,"\nEstimated Mu is: ",abs(res.const.bd$mu_par),"\n")
    
    
    if(Ntip(myTree)<5){
      cat("\nToo few tips to fit DDD and Time-varying models.\n")
      res.all <- as.data.frame(t(c("Clade.age"=tot_time,
                                   "CB.AICc"=res.const.b$aicc,
                                   "CB.lamb_par"=res.const.b$lamb_par,
                                   "CBD.AICc"=res.const.bd$aicc,
                                   "CBD.lamb_par"=res.const.bd$lamb_par,
                                   "CBD.mu_par"=res.const.bd$mu_par,
                                   "Best.DDD.model"=NA,
                                   "DDD.aicc"=NA,
                                   "DDD.lambda"=NA,
                                   "DDD.mu"=NA,
                                   "DDD.K"=NA,
                                   "Best Time Model"=NA,
                                   "Time.aicc"=NA,
                                   "Time.lambda"=NA,
                                   "Time.alpha"=NA,
                                   "Time.omega"=NA,
                                   "Time.mu"=NA,
                                   "Time.beta"=NA,
                                   "Time.pi"=NA)))
      cat("\nHence, LH values and AICc are calculated only for the constant-rate models.\n")
      return(res.all)
    }
    
    
    ##################### Peaks and dips in diversification ###########################
    ###################################################################################
    cat('\nEstimating positions for local peaks (if any)...\n')
    
    ###################### Function 1: Width of peak (or dip) constant ################
    cat("\nEstimating positions of peaks or dips keeping width of peak constant...\n")
    birth.epi.fix <- estimate_peak_or_dip_position(Tr,tot_time,res.const.bd,fix.width = TRUE,peak.in.which.rate = "birth")
    y.b.fix <- as.vector(unlist(c(birth.epi.fix[1,2:3])))
    p.b.fix <- as.vector(unlist(c(birth.epi.fix[1,1])))
    
    death.epi.fix <- estimate_peak_or_dip_position(Tr,tot_time,res.const.bd,fix.width = TRUE,peak.in.which.rate = "death")
    y.d.fix <- as.vector(unlist(c(death.epi.fix[1,2:3])))
    p.d.fix <- as.vector(unlist(c(death.epi.fix[1,1])))
    
    ############### Function 2: width of peak (or dip) to be optimized ################
    cat("\nEstimating positions of peaks or dips and letting width of peak vary...\n")
    birth.epi.var <- estimate_peak_or_dip_position(Tr,tot_time,res.const.bd,fix.width = FALSE,peak.in.which.rate = "birth")
    y.b.var <- as.vector(unlist(c(birth.epi.var[1,2:4])))
    p.b.var <- as.vector(unlist(c(birth.epi.var[1,1])))
    
    death.epi.var <- estimate_peak_or_dip_position(Tr,tot_time,res.const.bd,fix.width = FALSE,peak.in.which.rate = "death")
    y.d.var <- as.vector(unlist(c(death.epi.var[1,2:4])))
    p.d.var <- as.vector(unlist(c(death.epi.var[1,1])))
    
    
    ################################## Diversity-dependent models ####################################
    ##################################################################################################
    K <- Ntip(Tr) # expected carrying capacity
    
    lambda_ini <- abs(res.const.bd$lamb_par)
    mu_ini <- abs(res.const.bd$mu_par)
    
    cat("\nInitial parameters set for DDD models.\nCalculating likelihoods for DDD models......")
    
    DDD.models <- c("DDL + E","DDX + E","DD + EL","DD + EX")  # DDL + E worked for all the lineages. However, the remaining 3 did not converge in few cases
    #  in all the successful cases DDL + E had the best support
    LH.DDD <- df <- conv <- lamb.DDD <- mu.DDD <- K.DDD <- rep(NA,length(DDD.models))
    for(ddmodel in 1:length(DDD.models)){
      #ddmodel <- 1
      tryCatch(
        expr={
          res.DDD <- dd_ML(branching.times(Tr),c(lambda_ini,mu_ini,K),optimmethod = "simplex",ddmodel = ddmodel)
          cat("\nfitted this DDD model\n")
          LH.DDD[ddmodel] <- res.DDD$loglik
          df[ddmodel] <- res.DDD$df
          conv[ddmodel] <- res.DDD$conv
          lamb.DDD[ddmodel] <- res.DDD$lambda
          mu.DDD[ddmodel] <- res.DDD$mu
          K.DDD[ddmodel] <- res.DDD$K
        },
        error=function(e){},
        warning=function(w){},
        finally ={}
      )
      cat("\nDone with ",DDD.models[ddmodel],"\n")
    }
    cat("\nDDD LH estimates are:",LH.DDD,"\ndf:",df,"\nconv:",conv,"\n\n")
    
    
    l.conv <- length(which(!(conv %in% c(-1,NA))))
    if(l.conv < length(DDD.models)){
      df <- df[-which(conv == -1)]
      DDD.models <- DDD.models[-which(conv %in% c(-1,NA))]
      LH.DDD <- LH.DDD[-which(conv %in% c(-1,NA))]
      lamb.DDD <- lamb.DDD[-which(conv %in% c(-1,NA))]
      mu.DDD <- mu.DDD[-which(conv %in% c(-1,NA))]
      K.DDD <- K.DDD[-which(conv %in% c(-1,NA))]
    }
    
    if(length(df)==0){
      res.DDD <- as.data.frame(t(c("Best DDD model"=NA,"LH"=NA,"df"=NA,"aicc"=NA,
                                   "lambda"=NA,"mu"=NA,"K"=NA)))
    }else{
      if(Ntip(Tr)==(1+df+1)[1]){
        aicc.DDD <- 2*df - 2*LH.DDD
      }else{
        aicc.DDD <- 2*df - 2*LH.DDD + (2*df*(df+1))/(Ntip(Tr)-1-df-1)
      }
      
      names(aicc.DDD) <- DDD.models
      min.aicc.DDD <- min(aicc.DDD)
      max.df.DDD <- df[which(aicc.DDD == min.aicc.DDD)]
      max.LH.DDD <- LH.DDD[which(aicc.DDD == min.aicc.DDD)]
      best.model.DDD <- names(which(aicc.DDD == min.aicc.DDD))
      
      max.lamb.DDD <- lamb.DDD[which(aicc.DDD == min.aicc.DDD)]
      max.mu.DDD <- mu.DDD[which(aicc.DDD == min.aicc.DDD)]
      max.K.DDD <- K.DDD[which(aicc.DDD == min.aicc.DDD)]
      
      res.DDD <- as.data.frame(t(c("Best DDD model"=best.model.DDD,"LH"=max.LH.DDD,"df"=max.df.DDD ,"aicc"=min.aicc.DDD,
                                   "lambda"=max.lamb.DDD,"mu"=max.mu.DDD,"K"=max.K.DDD)))
    }
    
    
    
    
    ##################################################################################################
    #################################### Time-dependent models #######################################
    # THE NAMES OF THE MODELS HAVE BEEN CHANGED HERE FOR IN-CODE CHECKING WHILE MODEL FITTING IN THE COMING BIT OF THE CODE
    # THERE IS AN IF ELSE CONSTRUCT THAT UTILISES THE NAMES. MY NOMENCLATURE OF THE FUNCTIONS WAS EASIER TO WORK WITH 
    
    # Functions for speciation
    f.time.b.lin.abs <- function(t,y){y[1] + y[2] * t} # Equivalent to the BVARLIN function with absolute values
    f.time.b.lin.max <- function(t,y){sapply((y[1] + y[2] * t),max,0)} # Equivalent to the BVARLIN function with the maximum between the non-absolute values and zero. function suggested by morlon 2020
    f.time.b.exp.abs <- function(t,y){y[1] * exp(y[2] * t)} # Equivalent to the BVAREXP function with absolute values
    f.time.b.exp.max <- function(t,y){sapply((y[1] * exp(y[2] * t)),max,0)} # Equivalent to the BVAREXP function with the maximum between the non-absolute values and zero. function suggested by morlon 2020
    f.time.b.epi.w.fix <- function(t,y){y[1] + (y[2]*1/p.b.fix) * (t/((t-p.b.fix)^2+1))} # Function defined in this study to model symmetrical episodic rate shifts (birth) with a fixed width of width parameter = 1
    f.time.b.epi.w.var <- function(t,y){y[1] + (y[2]*abs(y[3])/p.b.var) * (t/((t-p.b.var)^2+abs(y[3])))} # Function defined in this study to model symmetrical episodic rate shifts (birth) with a varying width
    f.time.b.cst <- function(t,y){abs(y[1])} # Equivalent to the BCST function
    
    # Fucntions for extinctions
    f.time.d.lin.abs <- function(t,y){y[1] + y[2] * t} # Equivalent to the DVARLIN function with absolute values
    f.time.d.lin.max <- function(t,y){sapply((y[1] + y[2] * t),max,0)} # Equivalent to the DVARLIN function with the maximum between the non-absolute values and zero. function suggested by morlon 2020
    f.time.d.exp.abs <- function(t,y){y[1] * exp(y[2] * t)} # Equivalent to the DVAREXP function with absolute values
    f.time.d.exp.max <- function(t,y){sapply((y[1] * exp(y[2] * t)),max,0)} # Equivalent to the DVAREXP function with the maximum between the non-absolute values and zero. function suggested by morlon 2020
    f.time.d.epi.w.fix <- function(t,y){y[1] + (y[2]*1/p.d.fix) * (t/((t-p.d.fix)^2+1))} # Function defined in this study to model symmetrical episodic rate shifts (death) with a fixed width of width parameter = 1
    f.time.d.epi.w.var <- function(t,y){y[1] + (y[2]*abs(y[3])/p.d.var) * (t/((t-p.d.var)^2+abs(y[3])))} # Function defined in this study to model symmetrical episodic rate shifts (death) with a varying width
    f.time.d.cst <- function(t,y){abs(y[1])} # Equivalent to the DCST function
    f.time.d.0 <- function(t,y){0} # Equivalent to the D0 function
    
    # MAKING A NAMED VECTOR OF THE FUNCTIONS - FOR BIRTH RATES
    vector.b.funcs.time <- c("b.cst"=f.time.b.cst,
                             "b.exp.abs"=f.time.b.exp.abs,
                             "b.exp.max"=f.time.b.exp.max,
                             "b.lin.abs"=f.time.b.lin.abs,
                             "b.lin.max"=f.time.b.lin.max,
                             "b.epi.fix"=f.time.b.epi.w.fix,
                             "b.epi.var"=f.time.b.epi.w.var)
    names(vector.b.funcs.time)[c(6,7)] <- paste0(names(vector.b.funcs.time)[c(6,7)],".",c(p.b.fix,p.b.var))
    
    # MAKING A NAMED VECTOR OF THE FUNCTIONS - FOR DEATH RATES
    vector.d.funcs.time <- c("d.0"=f.time.d.0,
                             "d.cst"=f.time.d.cst,
                             "d.exp.abs"=f.time.d.exp.abs,
                             "d.exp.max"=f.time.d.exp.max,
                             "d.lin.abs"=f.time.d.lin.abs,
                             "d.lin.max"=f.time.d.lin.max,
                             "d.epi.fix"=f.time.d.epi.w.fix,
                             "d.epi.var"=f.time.d.epi.w.var)
    names(vector.d.funcs.time)[c(7,8)] <- paste0(names(vector.d.funcs.time)[c(7,8)],".",c(p.d.fix,p.d.var))
    
    
    
    cat("\nSetting initial parameters for Time-dependent model.",
        "\nCalculating likelihood for Time-dependent model......")
    
    
    LH.time.diff.funcs <- aicc.time.diff.funcs <- lamb.time <- alpha.time <- omega.time <- 
      mu.time <- beta.time <- pi.time <- complexity <- rep(NA,length(vector.b.funcs.time)*length(vector.d.funcs.time))
    res.count <- 1
    models.labels.time <- c()
    
    # running over all possible Time-dependent models for "Tr" phylogenetic tree
    for(m in 1:length(vector.b.funcs.time)){ # going over all b.funcs - speciation functions
      #m=2
      
      # setting initial parameters for lambda and algorithm
      if(substr(labels(vector.b.funcs.time)[m],3,5) == "cst"){
        lamb_par_init.time <- res.const.bd$lamb_par
        cst.lamb <- TRUE
        expo.lamb <- FALSE
      }else{
        lamb_par_init.time <- c(res.const.bd$lamb_par, runif(1,0,fin.lamb/2))
        
        cst.lamb <- FALSE
        if(substr(labels(vector.b.funcs.time)[m],3,5) == "exp"){
          expo.lamb <- TRUE
        }else{
          expo.lamb <- FALSE
        }
        
      }
      
      if(substr(labels(vector.b.funcs.time)[m],3,5) == "epi"){
        if(substr(labels(vector.b.funcs.time)[m],7,9) == "fix"){
          lamb_par_init.time <- y.b.fix
          if(is.na(p.b.fix)) next
        }else{
          lamb_par_init.time <- y.b.var
          if(is.na(p.b.var)) next
        }
        
      }
      
      # Time-dependent function for lambda in this iteration
      lamb.f <- vector.b.funcs.time[[m]]
      
      for(n in 1:length(vector.d.funcs.time)){ # going over all d.funcs - extinction functions
        #n=1
        if(m==1){
          if(n %in% 1:2) next
        }
        # setting initial parameters for mu and algorithm
        if(substr(labels(vector.d.funcs.time)[n],3,3) == "0"){
          mu_par_init.time <- c()
          fix.mu <- TRUE
          cst.mu <- TRUE
          expo.mu <- FALSE
        }else if(substr(labels(vector.d.funcs.time)[n],3,5) == "cst"){
          mu_par_init.time <- res.const.bd$mu_par
          cst.mu <- TRUE
          fix.mu <- FALSE
          expo.mu <- FALSE
        }else{
          mu_par_init.time <- c(res.const.bd$mu_par, runif(1,0,0.01))
          
          cst.mu <- FALSE
          fix.mu <- FALSE
          if(substr(labels(vector.d.funcs.time)[n],3,5) == "exp"){
            expo.mu <- TRUE
          }else{
            expo.mu <- FALSE
          }
        }
        
        if(substr(labels(vector.d.funcs.time)[n],3,5) == "epi"){
          if(substr(labels(vector.d.funcs.time)[n],7,9) == "fix"){
            mu_par_init.time <- y.d.fix
            if(is.na(p.d.fix)) next
          }else{
            mu_par_init.time <- y.d.var
            if(is.na(p.d.var)) next
          }
          
        }
        
        # Time-dependent function for mu in this iteration
        mu.f <- vector.d.funcs.time[[n]]
        
        error_counter <- 0      
        tryCatch(
          expr = {
            res.time <- fit_bd(Tr,max(branching.times(Tr)),lamb.f,mu.f,lamb_par_init.time,mu_par_init.time,f=1,dt=1e-3,
                               cst.lamb = cst.lamb,cst.mu = cst.mu,fix.mu = fix.mu, expo.lamb = expo.lamb,expo.mu = expo.mu)
            
            complexity[res.count] <- length(res.time$lamb_par) + length(res.time$mu_par)
            lamb.time[res.count] <- res.time$lamb_par[1]
            if(substr(labels(vector.b.funcs.time)[m],3,5) != "cst"){
              alpha.time[res.count] <- res.time$lamb_par[2]
            }
            if(substr(labels(vector.d.funcs.time)[n],3,3) != "0"){
              mu.time[res.count] <- res.time$mu_par[1]
              
              if(substr(labels(vector.d.funcs.time)[n],3,5) != "cst"){
                beta.time[res.count] <- res.time$mu_par[2]
              }
            }
            if(substr(labels(vector.b.funcs.time)[m],3,5) == "epi"){ ##### CHECK
              omega.time[res.count] <- res.time$lamb_par[3]
            }
            if(substr(labels(vector.d.funcs.time)[m],3,5) == "epi"){
              pi.time[res.count] <- res.time$mu_par[3]
            }
            
            LH.time.diff.funcs[res.count] <- res.time$LH
            
            if(res.time$aicc < 0){
              aicc.time.diff.funcs[res.count] <- NA
            }else{
              aicc.time.diff.funcs[res.count] <- res.time$aicc
            }
            
            
          },
          error = function(e){
            #cat("\nError   ",error_counter,"\n")
          },
          warning = function(w){},
          finally = {}
        )
        
        models.labels.time <- c(models.labels.time, paste0(labels(vector.b.funcs.time)[m],"_",labels(vector.d.funcs.time)[n]))
        cat("\nLikelihood for ",models.labels.time[res.count],", of complexity - ",complexity[res.count],", is:",LH.time.diff.funcs[res.count]," and AICc is:",aicc.time.diff.funcs[res.count],"\n")
        res.count = res.count + 1
      }
    }
    cat("\n\n\nRan all Time dependent models\n\n\n")
    
    
    # Choosing best model based on aicc scores and model complexity
    if(length(which(!is.na(aicc.time.diff.funcs))) != 0){
      
      min.aicc.time <- min(aicc.time.diff.funcs[!is.na(aicc.time.diff.funcs)]) # get the minimum AICc score
      cat("\nLeast AICc score:",min.aicc.time,"\n")
      
      # check if there are other models that aren't significantly worse than the one with lowest AICc
      ind.equi.models <- which(aicc.time.diff.funcs <= (min.aicc.time + 2)) # indices of equivalent models
      ind.best.model <- ind.equi.models
      
      eq.models.complexity <- complexity[ind.equi.models]
      eq.models.labels <- models.labels.time[ind.equi.models]
      eq.models.aicc <- aicc.time.diff.funcs[ind.equi.models]
      
      if(length(eq.models.labels)>1){
        cat("\nThere are multiple models with similar AICc scores.\n
          Complexity of best equivalent models -",eq.models.labels,":",eq.models.complexity,"respectively\n\n")
        
        ind.least.complexity <- which(models.labels.time %in% eq.models.labels[which(eq.models.complexity == min(eq.models.complexity))])
        ind.best.model <- ind.least.complexity
        
        # if complexities are also same then this section gets implemented
        if(length(ind.least.complexity)>1){
          #check if AICc values also exactly match
          ind.least.aicc <- which(models.labels.time %in% eq.models.labels[which(eq.models.aicc == min.aicc.time)])
          ind.best.model <- ind.least.aicc
          
          if(length(ind.least.aicc)>1){ # if AICc values of multiple models match exactly, then randomly choose one of the models, as all of them are equally likely
            ind.least.aicc <- ind.least.complexity[as.integer(runif(n = 1,min = 1,max = length(ind.least.complexity)))]
            ind.best.model <- ind.least.aicc
          }
        }
      }
      
      min.aicc.time <- aicc.time.diff.funcs[ind.best.model]
      max.LH.time <- LH.time.diff.funcs[ind.best.model]
      max.lamb.time <- lamb.time[ind.best.model]
      max.alpha.time <- alpha.time[ind.best.model]
      max.omega.time <- omega.time[ind.best.model]
      max.mu.time <- mu.time[ind.best.model]
      max.beta.time <- beta.time[ind.best.model]
      max.pi.time <- pi.time[ind.best.model]
      best.model.time <- models.labels.time[ind.best.model]
      
      cat("\nBest-fitting model:",best.model.time,"\n")
      res.time <- as.data.frame(t(c(best.model.time,max.lamb.time,max.alpha.time,max.omega.time,max.mu.time,max.beta.time,max.pi.time,min.aicc.time,max.LH.time)))
      colnames(res.time) <- c("Best time Model","lambda","alpha","omega","mu","beta","pi","aicc","LH")
    }else{
      res.time <- as.data.frame(t(c("Best time Model"=NA,"lambda"=NA,"alpha"=NA,"omega"=NA,"mu"=NA,"beta"=NA,"pi"=NA,"aicc"=NA,"LH"=NA)))
    }
    
    res.all <- as.data.frame(t(c("Clade.age"=tot_time,
                                 "CB.AICc"=res.const.b$aicc,
                                 "CB.lamb_par"=res.const.b$lamb_par,
                                 "CBD.AICc"=res.const.bd$aicc,
                                 "CBD.lamb_par"=res.const.bd$lamb_par,
                                 "CBD.mu_par"=res.const.bd$mu_par,
                                 "Best.DDD.model"=res.DDD$`Best DDD model`,
                                 "DDD.aicc"=res.DDD$aicc,
                                 "DDD.lambda"=res.DDD$lambda,
                                 "DDD.mu"=res.DDD$mu,
                                 "DDD.K"=res.DDD$K,
                                 "Best Time Model"=res.time$`Best time Model`,
                                 "Time.aicc"=res.time$aicc,
                                 "Time.lambda"=res.time$lambda,
                                 "Time.alpha"=res.time$alpha,
                                 "Time.omega"=res.time$omega,
                                 "Time.mu"=res.time$mu,
                                 "Time.beta"=res.time$beta,
                                 "Time.pi"=res.time$pi)))
    cat("\nCalculated LH values and AICc\n")
    return(res.all)
  }
  
  
  # Function to extract largest sub-clade which is younger than the maximum age of a certain paleo-climatic variable
  extract_largest_subclades_numbers <- function(env_data,myTree){
    max.Env.Age <- max(env_data$Age)
    nodes_ages <- tree.age(myTree)
    
    cr.age <- max(nodes_ages$ages)
    
    if(cr.age > max.Env.Age){
      # extracting node numbers and ages of all nodes
      
      nodes_ages <- nodes_ages[-which(nodes_ages$ages == 0.000),]
      nodes_ages <- nodes_ages[-which(nodes_ages$ages > max.Env.Age),]
      if(length(nodes_ages[,1]) == 0){
        cat("\nThere aren't any clade younger than the maximum age of the paleo-climate data\n")
        return(NA)
      }
      clades.node.nums <- c()
      
      while(TRUE){
        mx <- as.numeric(nodes_ages$elements[which(nodes_ages$ages == max(nodes_ages$ages))])
        d.mx <- Descendants(myTree,mx,"all")
        intersecting.d.mx <- intersect(as.numeric(nodes_ages$elements),d.mx)
        if(length(intersecting.d.mx) != 0){
          for(des in 1:length(intersecting.d.mx)){
            nodes_ages <- nodes_ages[-which(nodes_ages$elements == intersecting.d.mx[des]),]
          }
        }
        nodes_ages <- nodes_ages[-which(nodes_ages$elements == mx),]
        clades.node.nums <- c(clades.node.nums,mx)
        if(length(nodes_ages[,1]) == 0){break}
        else{next}
      }
    }else{
      clades.node.nums <- as.numeric(c(nodes_ages$elements[which(nodes_ages$ages == max(nodes_ages$ages))]))
    }
    return(clades.node.nums)
  }
  
  
  model.comparisons.wo.env <- matrix(ncol = 21)
  colnames(model.comparisons.wo.env)<- c("lineage.names","N.tips","Clade.age","CB.AICc","CB.lamb_par","CBD.AICc","CBD.lamb_par","CBD.mu_par",
                                         "Best.DDD.model","DDD.aicc","DDD.lambda","DDD.mu","DDD.K",
                                         "Best Time Model","Time.aicc","Time.lambda","Time.alpha","Time.omega","Time.mu","Time.beta","Time.pi")
  
  # The following matrix is for storing parameter estimates when we just consider temperature as the only environmental variable to compare diversity-dependence with
  model.comparisons.w.temp <- matrix(ncol = 28)
  colnames(model.comparisons.w.temp) <- c("lineage.names","clades.node.nums","N.tips","Clade.age","CB.AICc","CB.lamb_par","CBD.AICc","CBD.lamb_par","CBD.mu_par",
                                          "Best.DDD.model","DDD.aicc","DDD.lambda","DDD.mu","DDD.K",
                                          "Best Time Model","Time.aicc","Time.lambda","Time.alpha","Time.omega","Time.mu","Time.beta","Time.pi",
                                          "Best Temp Model","Temp.aicc","Temp.lambda","Temp.alpha","Temp.mu","Temp.beta")
  
  model.comparisons.w.him.temp <- matrix(ncol = 34)
  colnames(model.comparisons.w.him.temp) <- c("lineage.names","clades.node.nums","N.tips","Clade.age","CB.AICc","CB.lamb_par","CBD.AICc","CBD.lamb_par","CBD.mu_par",
                                              "Best.DDD.model","DDD.aicc","DDD.lambda","DDD.mu","DDD.K",
                                              "Best Time Model","Time.aicc","Time.lambda","Time.alpha","Time.omega","Time.mu","Time.beta","Time.pi",
                                              "Best Temp Model","Temp.aicc","Temp.lambda","Temp.alpha","Temp.mu","Temp.beta",
                                              "Best Him.oro Model","Him.oro.aicc","Him.oro.lambda","Him.oro.alpha","Him.oro.mu","Him.oro.beta")
  
  
  model.comparisons.w.c4.him.temp <- matrix(ncol = 40)
  colnames(model.comparisons.w.c4.him.temp) <- c("lineage.names","clades.node.nums","N.tips","Clade.age","CB.AICc","CB.lamb_par","CBD.AICc","CBD.lamb_par","CBD.mu_par",
                                                 "Best.DDD.model","DDD.aicc","DDD.lambda","DDD.mu","DDD.K",
                                                 "Best Time Model","Time.aicc","Time.lambda","Time.alpha","Time.omega","Time.mu","Time.beta","Time.pi",
                                                 "Best Temp Model","Temp.aicc","Temp.lambda","Temp.alpha","Temp.mu","Temp.beta",
                                                 "Best Him.oro Model","Him.oro.aicc","Him.oro.lambda","Him.oro.alpha","Him.oro.mu","Him.oro.beta",
                                                 "Best C4-Exp Model","C4-Exp.aicc","C4-Exp.lambda","C4-Exp.alpha","C4-Exp.mu","C4-Exp.beta")
  #lineage.name <- lineage.names[3]
  
  myTree <- read.tree(paste0(lineage.name,".tre"))
  #system2("mkdir",paste0("../model_comparisons_parallelized/bc_",lineage.name))
  
  
  cr.age <- max(branching.times(myTree))
  cat(paste0("\nCurrently using tree of ",lineage.name,"\n"))
  
  ################### NON-PALEOCLIMATIC_VARIABLE DEPENDENT MODELS ################
  # get the best paramters for CB, CBD, DDD and tD (time-dependence)
  res.npc1 <- run_model_fit_wo_env(myTree)
  
  # keep storing those values in this matrix 
  model.comparisons.wo.env <- rbind(model.comparisons.wo.env,
                                    c(lineage.name,
                                      Ntip(myTree),
                                      res.npc1$Clade.age,
                                      res.npc1$CB.AICc,
                                      res.npc1$CB.lamb_par,
                                      res.npc1$CBD.AICc,
                                      res.npc1$CBD.lamb_par,
                                      res.npc1$CBD.mu_par,
                                      res.npc1$Best.DDD.model,
                                      res.npc1$DDD.aicc,
                                      res.npc1$DDD.lambda,
                                      res.npc1$DDD.mu,
                                      res.npc1$DDD.K,
                                      res.npc1$`Best Time Model`,
                                      res.npc1$Time.aicc,
                                      res.npc1$Time.lambda,
                                      res.npc1$Time.alpha,
                                      res.npc1$Time.omega,
                                      res.npc1$Time.mu,
                                      res.npc1$Time.beta,
                                      res.npc1$Time.pi))
  
  # since it every lineage is parallely analysed, we are saving the result for each lineage as a separate file, 
  # all of these files would be clubbed later
  # writing the result of fitting models without environmental dependence
  write.csv(model.comparisons.wo.env[-1,],file = paste0("../model_comparisons_parallelized/model_comparisons_w0_env_model_concatenated_very_new_",lineage.name,".csv"))
  ############################################################################################
  
  
  
  ############################## Temperature-dependent models #########################
  # Functions for speciation
  f.temp.b.lin.abs <- function(t,x,y){y[1] + y[2] * x}
  f.temp.b.lin.max <- function(t,x,y){sapply((y[1] + y[2] * x),max,0)}
  f.temp.b.exp.abs <- function(t,x,y){y[1] * exp(y[2] * x)}
  f.temp.b.exp.max <- function(t,x,y){sapply(y[1] * exp(y[2] * x),max,0)}
  f.temp.b.cst <- function(t,x,y){y[1]}
  f.temp.b.mtb.abs <- function(t,x,y){y[1] * exp(-y[2] / x)}
  f.temp.b.mtb.max <- function(t,x,y){sapply(y[1] * exp(-y[2] / x),max,0)}
  
  # Fucntions for extinctions
  f.temp.d.lin.abs <- function(t,x,y){y[1] + y[2] * x}
  f.temp.d.lin.max <- function(t,x,y){sapply((y[1] + y[2] * x),max,0)}
  f.temp.d.exp.abs <- function(t,x,y){y[1] * exp(y[2] * x)}
  f.temp.d.exp.max <- function(t,x,y){sapply(y[1] * exp(y[2] * x),max,0)}
  f.temp.d.cst <- function(t,x,y){y[1]}
  f.temp.d.0 <- function(t,x,y){0}
  f.temp.d.mtb.abs <- function(t,x,y){y[1] * exp(-y[2] / x)}
  f.temp.d.mtb.max <- function(t,x,y){sapply(y[1] * exp(-y[2] / x),max,0)}
  
  vector.b.funcs.temp <- c("b.cst"=f.temp.b.cst,
                           "b.exp.abs"=f.temp.b.exp.abs,
                           "b.exp.max"=f.temp.b.exp.max,
                           "b.lin.abs"=f.temp.b.lin.abs,
                           "b.lin.max"=f.temp.b.lin.max,
                           "b.mtb.abs"=f.temp.b.mtb.abs,
                           "b.mtb.max"=f.temp.b.mtb.max)
  
  vector.d.funcs.temp <- c("d.0"=f.temp.d.0,
                           "d.cst"=f.temp.d.cst,
                           "d.exp.abs"=f.temp.d.exp.abs,
                           "d.exp.max"=f.temp.d.exp.max,
                           "d.lin.abs"=f.temp.d.lin.abs,
                           "d.lin.max"=f.temp.d.lin.max,
                           "d.mtb.abs"=f.temp.d.mtb.abs,
                           "d.mtb.max"=f.temp.d.mtb.max)
  
  ############################## Himalayan-Orogeny-dependent models #################################
  ##################################################################################################
  # Functions for speciation
  f.him.oro.b.lin.abs <- function(t,x,y){y[1] + y[2] * x}
  f.him.oro.b.lin.max <- function(t,x,y){sapply((y[1] + y[2] * x),max,0)}
  f.him.oro.b.exp.abs <- function(t,x,y){y[1] * exp(y[2] * x)}
  f.him.oro.b.exp.max <- function(t,x,y){sapply(y[1] * exp(y[2] * x),max,0)}
  f.him.oro.b.cst <- function(t,x,y){y[1]}
  
  # Fucntions for extinctions
  f.him.oro.d.lin.abs <- function(t,x,y){y[1] + y[2] * x}
  f.him.oro.d.lin.max <- function(t,x,y){sapply((y[1] + y[2] * x),max,0)}
  f.him.oro.d.exp.abs <- function(t,x,y){y[1] * exp(y[2] * x)}
  f.him.oro.d.exp.max <- function(t,x,y){sapply(y[1] * exp(y[2] * x),max,0)}
  f.him.oro.d.cst <- function(t,x,y){y[1]}
  f.him.oro.d.0 <- function(t,x,y){0}
  
  vector.b.funcs.him.oro <- c("b.cst"=f.him.oro.b.cst,
                              "b.exp.abs"=f.him.oro.b.exp.abs,
                              "b.exp.max"=f.him.oro.b.exp.max,
                              "b.lin.abs"=f.him.oro.b.lin.abs,
                              "b.lin.max"=f.him.oro.b.lin.max)
  
  vector.d.funcs.him.oro <- c("d.0"=f.him.oro.d.0,
                              "d.cst"=f.him.oro.d.cst,
                              "d.exp.abs"=f.him.oro.d.exp.abs,
                              "d.exp.max"=f.him.oro.d.exp.max,
                              "d.lin.abs"=f.him.oro.d.lin.abs,
                              "d.lin.max"=f.him.oro.d.lin.max)
  
  
  # ############################# C4-Plant-Expansion-dependent models ################################
  # ##################################################################################################
  # Functions for speciation
  f.c4.exp.b.lin.abs <- function(t,x,y){y[1] + y[2] * x}
  f.c4.exp.b.lin.max <- function(t,x,y){sapply((y[1] + y[2] * x),max,0)}
  f.c4.exp.b.exp.abs <- function(t,x,y){y[1] * exp(y[2] * x)}
  f.c4.exp.b.exp.max <- function(t,x,y){sapply(y[1] * exp(y[2] * x),max,0)}
  f.c4.exp.b.cst <- function(t,x,y){y[1]}
  
  # Fucntions for extinctions
  f.c4.exp.d.lin.abs <- function(t,x,y){y[1] + y[2] * x}
  f.c4.exp.d.lin.max <- function(t,x,y){sapply((y[1] + y[2] * x),max,0)}
  f.c4.exp.d.exp.abs <- function(t,x,y){y[1] * exp(y[2] * x)}
  f.c4.exp.d.exp.max <- function(t,x,y){sapply(y[1] * exp(y[2] * x),max,0)}
  f.c4.exp.d.cst <- function(t,x,y){y[1]}
  f.c4.exp.d.0 <- function(t,x,y){0}
  
  vector.b.funcs.c4.exp <- c("b.cst"=f.c4.exp.b.cst,
                             "b.exp.abs"=f.c4.exp.b.exp.abs,
                             "b.exp.max"=f.c4.exp.b.exp.max,
                             "b.lin.abs"=f.c4.exp.b.lin.abs,
                             "b.lin.max"=f.c4.exp.b.lin.max)
  
  vector.d.funcs.c4.exp <- c("d.0"=f.c4.exp.d.0,
                             "d.cst"=f.c4.exp.d.cst,
                             "d.exp.abs"=f.c4.exp.d.exp.abs,
                             "d.exp.max"=f.c4.exp.d.exp.max,
                             "d.lin.abs"=f.c4.exp.d.lin.abs,
                             "d.lin.max"=f.c4.exp.d.lin.max)
  
  
  # Don't continue to check for environmental variable dependence if the number of tips in this lineage is less than 10
  if(Ntip(myTree)<10){
    cat("\n\nDone with ",lineage.name,"\n\n")
    return()
  }
  
  # In this following section, plaeo-climate-data-dependence is checked.
  # The maximum age till the paleodata is available is a major constraint as it won't fit a lineage that is older than the maximium age in the paleodata.
  # Hence we need to break our phylogeny of interest into sub-clades and take the largest sub-clades that are younger than the maximum age in the paleodata.
  # Also, we are considering only those clades which have more than 10 tips for checking the influence of environmental variables in their diversification.
  # We chose 10 tips as a filter because in our pilot analyses clades with less than 10 tips were showing errors while fitting the environmental models.
  
  
  # However for the final inference of the study we only used the clades as a whole, and not specific sub-clades within them to check the influence of
  # an environmental variable. We have however retained the code that involves breaking a clade of interest into sub-clades which are younger than the
  # age of an environmental variable
  ##################################### ONLY TEMPERATURE DEPENDENCE ###################################
  cat("\nChecking temperature dependence.....\n")
  
  clades.node.nums <- extract_largest_subclades_numbers(Temp_dat,myTree)
  if(is.na(clades.node.nums)[1]) {cat("Done with ",lineage.name,"\n\n");return()}  # I am confidently putting return() here as the Temperature-data has the oldest date of all
  # so ignoring himalayan-orogent and c4-expansion is okay
  
  cat("\n\nClades no.s to be checked:",clades.node.nums,"\n\n")
  
  # Going through all largest subtrees
  for(tr.i in 1:length(clades.node.nums)){
    #tr.i = 1
    
    Tr <- extract.clade(phy = myTree,node = clades.node.nums[tr.i])
    #Ntip(Tr)
    if(Ntip(Tr)<10){
      cat("\nNumber of tips is too less for this clade for performing analysis.\n\n")
      next
    }
    
    times <- branching.times(Tr)
    tot_time<-max(times)
    
    sorted.times <- rev(as.numeric(sort(times)))
    
    # Let's calculate the initial speciation rate (dummy) - rate at crown
    ini.lamb <- (1/(sorted.times[1]-sorted.times[2]))/3
    # # Let's calculate the final speciation rate (dummy) - rate at tips
    fin.lamb <- (1/(sorted.times[length(sorted.times)-1] - sorted.times[length(sorted.times)]))/Ntip(Tr)
    
    ############## Constant rate, Diversity Dependent, Time dependent models ##############
    res.npc2 <- run_model_fit_wo_env(Tr)
    res.const.bd <- as.data.frame(t(c("aicc"=res.npc2$CBD.AICc,"lamb_par"=res.npc2$CBD.lamb_par,"mu_par"=res.npc2$CBD.mu_par)))
    
    
    ############# temperature-dependence
    
    res.temp <- run_fit_env(Temp_dat,Tr,lineage.name,"Temp",vector.b.funcs.temp,vector.d.funcs.temp,res.const.bd)
    
    model.comparisons.w.temp <- rbind(model.comparisons.w.temp,
                                      c(lineage.name,
                                        clades.node.nums[tr.i],
                                        Ntip(Tr),
                                        res.npc2$Clade.age,
                                        res.npc2$CB.AICc,
                                        res.npc2$CB.lamb_par,
                                        res.npc2$CBD.AICc,
                                        res.npc2$CBD.lamb_par,
                                        res.npc2$CBD.mu_par,
                                        res.npc2$Best.DDD.model,
                                        res.npc2$DDD.aicc,
                                        res.npc2$DDD.lambda,
                                        res.npc2$DDD.mu,
                                        res.npc2$DDD.K,
                                        res.npc2$`Best Time Model`,
                                        res.npc2$Time.aicc,
                                        res.npc2$Time.lambda,
                                        res.npc2$Time.alpha,
                                        res.npc2$Time.omega,
                                        res.npc2$Time.mu,
                                        res.npc2$Time.beta,
                                        res.npc2$Time.pi,
                                        res.temp$`Best Env model`,
                                        res.temp$aicc,
                                        res.temp$lambda,
                                        res.temp$alpha,
                                        res.temp$mu,
                                        res.temp$beta))
  }
  
  # writing the result of fitting models with time + diversity + temperature dependence
  write.csv(model.comparisons.w.temp[-1,],file = paste0("../model_comparisons_parallelized/model_comparisons_w_smoothened_temp_model_concatenated_new_",lineage.name,".csv"))
  
  ################################## TEMPERATURE AND HIMALAYAN OROGENY #######################################
  cat("\nComparing temperature and himalayan-orogeny dependence.....\n")
  
  clades.node.nums <- extract_largest_subclades_numbers(Him.oro_dat,myTree)
  if(is.na(clades.node.nums)[1]) {cat("Done with ",lineage.name,"\n\n");return()}  # I am confidently putting "return()" here as the Himalayan-orogeny-data is of older date than c4-data
  # so ignoring c4-expansion is okay
  
  cat("\n\nClades no.s to be checked:",clades.node.nums,"\n\n")
  
  # Going through all largest subtrees
  for(tr.i in 1:length(clades.node.nums)){
    #tr.i = 1
    
    Tr <- extract.clade(phy = myTree,node = clades.node.nums[tr.i])
    #Ntip(Tr)
    if(Ntip(Tr)<10){
      cat("\nNumber of tips is too less for this clade for performing analyses.\n\n")
      next
    }
    
    times <- branching.times(Tr)
    tot_time<-max(times)
    
    sorted.times <- rev(as.numeric(sort(times)))
    
    # Let's calculate the initial speciation rate (dummy) - rate at crown
    ini.lamb <- (1/(sorted.times[1]-sorted.times[2]))/3
    # # Let's calculate the final speciation rate (dummy) - rate at tips
    fin.lamb <- (1/(sorted.times[length(sorted.times)-1] - sorted.times[length(sorted.times)]))/Ntip(Tr)
    
    ############## Constant rate, Diversity Dependent, Time dependent models ##############
    res.npc2 <- run_model_fit_wo_env(Tr)
    res.const.bd <- as.data.frame(t(c("aicc"=res.npc2$CBD.AICc,"lamb_par"=res.npc2$CBD.lamb_par,"mu_par"=res.npc2$CBD.mu_par)))
    
    
    ############# temperature-dependence
    res.temp <- run_fit_env(Temp_dat,Tr,lineage.name,"Temp",vector.b.funcs.temp,vector.d.funcs.temp,res.const.bd)
    
    ############### himalayan-orogeny dependence
    res.him.oro <- run_fit_env(Him.oro_dat,Tr,lineage.name,"Him.oro",vector.b.funcs.him.oro,vector.d.funcs.him.oro,res.const.bd)
    
    model.comparisons.w.him.temp <- rbind(model.comparisons.w.him.temp,
                                          c(lineage.name,
                                            clades.node.nums[tr.i],
                                            Ntip(Tr),
                                            res.npc2$Clade.age,
                                            res.npc2$CB.AICc,
                                            res.npc2$CB.lamb_par,
                                            res.npc2$CBD.AICc,
                                            res.npc2$CBD.lamb_par,
                                            res.npc2$CBD.mu_par,
                                            res.npc2$Best.DDD.model,
                                            res.npc2$DDD.aicc,
                                            res.npc2$DDD.lambda,
                                            res.npc2$DDD.mu,
                                            res.npc2$DDD.K,
                                            res.npc2$`Best Time Model`,
                                            res.npc2$Time.aicc,
                                            res.npc2$Time.lambda,
                                            res.npc2$Time.alpha,
                                            res.npc2$Time.omega,
                                            res.npc2$Time.mu,
                                            res.npc2$Time.beta,
                                            res.npc2$Time.pi,
                                            res.temp$`Best Env model`,
                                            res.temp$aicc,
                                            res.temp$lambda,
                                            res.temp$alpha,
                                            res.temp$mu,
                                            res.temp$beta,
                                            res.him.oro$`Best Env model`,
                                            res.him.oro$aicc,
                                            res.him.oro$lambda,
                                            res.him.oro$alpha,
                                            res.him.oro$mu,
                                            res.him.oro$beta))
  }
  
  # writing the result of fitting models with time + diversity + temperature + himalayan-orogeny dependence
  write.csv(model.comparisons.w.him.temp[-1,],file = paste0("../model_comparisons_parallelized/model_comparisons_w_smoothened_him_temp_model_concatenated_very_new_",lineage.name,".csv"))
  
  ################################## TEMPERATURE, HIMALAYAN OROGENY AND C4-PLANT EXPANSION #######################################
  cat("\nComparing temperature, himalayan-orogeny and c4-plant expansion dependence.....\n")
  clades.node.nums <- extract_largest_subclades_numbers(C_dat_potwar,myTree)
  if(is.na(clades.node.nums)[1]) {cat("Done with ",lineage.name,"\n\n");return()}
  
  cat("\n\nClades no.s to be checked:",clades.node.nums,"\n\n")
  
  # Going through all largest subtrees
  for(tr.i in 1:length(clades.node.nums)){
    #tr.i = 1
    
    Tr <- extract.clade(phy = myTree,node = clades.node.nums[tr.i])
    #Ntip(Tr)
    if(Ntip(Tr)<10){
      cat("\nNumber of tips is too less for this clade for performing analysis.\n\n")
      next
    }
    
    times <- branching.times(Tr)
    tot_time<-max(times)
    
    sorted.times <- rev(as.numeric(sort(times)))
    
    # Let's calculate the initial speciation rate (dummy) - rate at crown
    ini.lamb <- (1/(sorted.times[1]-sorted.times[2]))/3
    # # Let's calculate the final speciation rate (dummy) - rate at tips
    fin.lamb <- (1/(sorted.times[length(sorted.times)-1] - sorted.times[length(sorted.times)]))/Ntip(Tr)
    
    ############## Constant rate, Diversity Dependent, Time dependent models ##############
    res.npc2 <- run_model_fit_wo_env(Tr)
    res.const.bd <- as.data.frame(t(c("aicc"=res.npc2$CBD.AICc,"lamb_par"=res.npc2$CBD.lamb_par,"mu_par"=res.npc2$CBD.mu_par)))
    
    
    ############# temperature-dependence
    res.temp <- run_fit_env(Temp_dat,Tr,lineage.name,"Temp",vector.b.funcs.temp,vector.d.funcs.temp,res.const.bd)
    
    ############### himalayan-orogeny dependence
    res.him.oro <- run_fit_env(Him.oro_dat,Tr,lineage.name,"Him.oro",vector.b.funcs.him.oro,vector.d.funcs.him.oro,res.const.bd)
    
    ############### c4-plant expansion dependence
    res.c4.exp <- run_fit_env(C_dat_potwar,Tr,lineage.name,"C4-Exp",vector.b.funcs.c4.exp,vector.d.funcs.c4.exp,res.const.bd)
    
    model.comparisons.w.c4.him.temp <- rbind(model.comparisons.w.c4.him.temp,
                                             c(lineage.name,
                                               clades.node.nums[tr.i],
                                               Ntip(Tr),
                                               res.npc2$Clade.age,
                                               res.npc2$CB.AICc,
                                               res.npc2$CB.lamb_par,
                                               res.npc2$CBD.AICc,
                                               res.npc2$CBD.lamb_par,
                                               res.npc2$CBD.mu_par,
                                               res.npc2$Best.DDD.model,
                                               res.npc2$DDD.aicc,
                                               res.npc2$DDD.lambda,
                                               res.npc2$DDD.mu,
                                               res.npc2$DDD.K,
                                               res.npc2$`Best Time Model`,
                                               res.npc2$Time.aicc,
                                               res.npc2$Time.lambda,
                                               res.npc2$Time.alpha,
                                               res.npc2$Time.omega,
                                               res.npc2$Time.mu,
                                               res.npc2$Time.beta,
                                               res.npc2$Time.pi,
                                               res.temp$`Best Env model`,
                                               res.temp$aicc,
                                               res.temp$lambda,
                                               res.temp$alpha,
                                               res.temp$mu,
                                               res.temp$beta,
                                               res.him.oro$`Best Env model`,
                                               res.him.oro$aicc,
                                               res.him.oro$lambda,
                                               res.him.oro$alpha,
                                               res.him.oro$mu,
                                               res.him.oro$beta,
                                               res.c4.exp$`Best Env model`,
                                               res.c4.exp$aicc,
                                               res.c4.exp$lambda,
                                               res.c4.exp$alpha,
                                               res.c4.exp$mu,
                                               res.c4.exp$beta))
    
  }
  
  # writing the result of fitting models with time + diversity + temperature + himalayan-orogeny + c4-plant expansion dependence
  write.csv(model.comparisons.w.c4.him.temp[-1,],file = paste0("../model_comparisons_parallelized/model_comparisons_w_smoothened_c4_him_temp_model_concatenated_very_new_",lineage.name,".csv"))
  ###################################### DONE WITH MODELS ###########################################
  
  
  cat("Done with ",lineage.name,"\n\n")
}

do <- function(pop, fun, ncores = 5, ...) {
  require(foreach)
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  foreach(i = pop) %dopar% fun(i, ...)
}
do(lineage.names,main_func)
main_func(lineage.names[3])
# parLapply(cl,lineage.names,main_func)
# stopCluster(cl)

out.files <- list.files("../model_comparisons_parallelized/")
ind.c4.him.temp <- grep("smoothened_c4",out.files)
ind.him.temp <- grep("smoothened_him",out.files)
ind.temp <- grep("smoothened_temp",out.files)
ind.wo.env <- grep("w0",out.files)

dfs.wo.env <- lapply(X = paste0("../model_comparisons_parallelized/",out.files[ind.wo.env]),read.csv)
dfs.temp <- lapply(X = paste0("../model_comparisons_parallelized/",out.files[ind.temp]),read.csv)
dfs.him.temp <- lapply(X = paste0("../model_comparisons_parallelized/",out.files[ind.him.temp]),read.csv)
dfs.c4.him.temp <- lapply(X = paste0("../model_comparisons_parallelized/",out.files[ind.c4.him.temp]),read.csv)

# function to combine the results of each distinct analysis
combine_data <- function(dfs){
  out.matrix <- matrix(ncol=nrow(dfs[[1]]))
  colnames(out.matrix) <- dfs[[1]][,1]
  
  for(i in 1:length(dfs)){
    if(length(dfs[[i]]) != 2) next
    out.matrix <- rbind(out.matrix,dfs[[i]][,2])
  }
  return(out.matrix)
}

# combining those files
combined.wo.env <- combine_data(dfs.wo.env)
combined.temp <- combine_data(dfs.temp)
combined.him <- combine_data(dfs.him.temp)
combined.c4 <- combine_data(dfs.c4.him.temp)

system2("mkdir","../model_comparisons_parallelized/temp_files")
system2("mv",c("../model_comparisons_parallelized/*.csv","../model_comparisons_parallelized/temp_files"))

# Write all the matrices
write.csv(combined.temp[-1,],file = "../model_comparisons_parallelized/model_comparisons_w_smoothened_temp_model_concatenated_new.csv")
write.csv(combined.him[-1,],file = "../model_comparisons_parallelized/model_comparisons_w_smoothened_him_temp_model_concatenated_very_new.csv")
write.csv(combined.c4[-1,],file = "../model_comparisons_parallelized/model_comparisons_w_smoothened_c4_him_temp_model_concatenated_very_new.csv")

write.csv(combined.wo.env[-1,],file = "../model_comparisons_parallelized/model_comparisons_w0_env_model_concatenated_very_new.csv")
#######################################################################################################



##############################################################################################################
####################################### AICc-BASED BEST-FITTING MODEL-SELECTION ##############################

# code for selecting best model among the entire pool
df <- read.csv("../model_comparisons_parallelized/model_comparisons_w0_env_model_concatenated_very_new.csv")
#View(df)

# order it based on lineage names in alphabetical order
df <- df[order(df$lineage.names,decreasing=FALSE),]
nms <- colnames(df) # names of columns of df

aicc.cols <- nms[c(which(grepl("AICc",nms)),which(grepl("aicc",nms)))] # extract just the columns that have AICc scores
aicc.cols <- aicc.cols[-which(aicc.cols == "DDD.aicc")]

df <- cbind(df,"lowest.aicc"=as.vector(apply(df[,aicc.cols], 1, min, na.rm = TRUE))) # add a column for the lowest AICc score per lineage
df <- cbind(df,"delta.aicc.CB"=(df[,aicc.cols[1]]-df[,"lowest.aicc"])) # add a column for delta AICc for CB model
df <- cbind(df,"delta.aicc.CBD"=(df[,aicc.cols[2]]-df[,"lowest.aicc"])) # add a column for delta AICc for CBD model
df <- cbind(df,"delta.aicc.Time"=(df[,aicc.cols[3]]-df[,"lowest.aicc"])) # # add a column for delta AICc for Time-varying model

best.model <- vector("character",nrow(df))
best.model.paramters <- matrix(nrow=1,ncol=6)
colnames(best.model.paramters) <- c("lambda0","alpha","omega","mu0","beta","pi")

# Run through each lineage and check what's the best model among - CB, CBD and Time-varying
for(i in 1:nrow(df)){
  models.w.eq.aicc <- aicc.cols[which(df[i,(ncol(df)-2):ncol(df)]<2)]
  if(aicc.cols[1] %in% models.w.eq.aicc){
    best.model[i] <- 'CB'
    best.model.paramters <- rbind(best.model.paramters,c(df$CB.lamb_par[i],0,NA,0,0,NA))
    next
  }else if(aicc.cols[2] %in% models.w.eq.aicc){
    best.model[i] <- 'CBD'
    best.model.paramters <- rbind(best.model.paramters,c(df$CBD.lamb_par[i],0,NA,df$CBD.mu_par[i],0,NA))
    next
  }else{
    best.model[i] <- 'Time'
    best.model.paramters <- rbind(best.model.paramters,c(df$Time.lambda[i],df$Time.alpha[i],df$Time.omega[i],df$Time.mu[i],df$Time.beta[i],df$Time.pi[i]))
  }
  
}

best.model.paramters <- best.model.paramters[-1,] # remove the first row of the matrix as it would just be "NA"s
best.model.paramters[which(is.na(best.model.paramters))] <- 0 # replace all the other "NA"s in your matrix with 0-s

df <- as.data.frame(cbind(df,"best.model"=best.model,best.model.paramters)) # add a column for the names of the best model out of : CB, CBD, Time
ind.time.best <- which(df$best.model=="Time") # get the indices for all the lineages/rows where "Time" is the best model
df[ind.time.best,"best.model"] <- df[ind.time.best,"Best.Time.Model"] # replace "Time" with the actual time-varying model name in this section of the matrix
#View(df)
write.csv(df,"../model_comparisons_parallelized/new_table_wo_env3.csv")
##############################################################################################################################



########################################################################################################################################
################################ CALCULATE RATES THROUGH TIME USING THE BEST MODEL PARAMETERS ##########################################
df <- read.csv("../model_comparisons_parallelized/new_table_wo_env3.csv")
clade.ages <- df$Clade.age

df <- as.data.frame(cbind(df$lineage.names,
                          df$best.model,
                          df$lambda0,
                          df$alpha,
                          df$omega,
                          df$mu0,
                          df$beta,
                          df$pi))

colnames(df) <- c("Lineages","Best.model","lambda0","alpha","omega","mu0","beta","pi")

# make matrix to store 100 time slices for each lineage
# here each column will represent each lineage
timeseries_for_rates <- matrix(ncol=nrow(df),nrow=100)
for(i in 1:length(clade.ages)){
  #i=1
  timeseries_for_rates[,i] <- seq(clade.ages[i],clade.ages[i]/100,-clade.ages[i]/100)
}
colnames(timeseries_for_rates) <- df$Lineages


df <- df[order(df$Lineages,decreasing=FALSE),]
timeseries_for_rates <- timeseries_for_rates[,order(colnames(timeseries_for_rates),decreasing = FALSE)]

# Functions for speciation
f.time.b.lin.abs <- function(t,y){abs(y[1] + y[2] * t)}
f.time.b.lin.max <- function(t,y){abs(sapply((y[1] + y[2] * t),max,0))}
f.time.b.exp.abs <- function(t,y){abs(y[1] * exp(y[2] * t))}
f.time.b.exp.max <- function(t,y){abs(sapply((y[1] * exp(y[2] * t)),max,0))}
f.time.b.epi.w.fix <- function(t,y){abs(y[1] + (y[2]*1/p.b) * (t/((t-p.b)^2+1)))}
f.time.b.epi.w.var <- function(t,y){abs(y[1] + (y[2]*abs(y[3])/p.b) * (t/((t-p.b)^2+abs(y[3]))))}
f.time.b.cst <- function(t,y){rep(abs(y[1]),length(t))}

# Fucntions for extinctions
f.time.d.lin.abs <- function(t,y){abs(y[1] + y[2] * t)}
f.time.d.lin.max <- function(t,y){abs(sapply((y[1] + y[2] * t),max,0))}
f.time.d.exp.abs <- function(t,y){abs(y[1] * exp(y[2] * t))}
f.time.d.exp.max <- function(t,y){abs(sapply((y[1] * exp(y[2] * t)),max,0))}
f.time.d.epi.w.fix <- function(t,y){abs(y[1] + (y[2]*1/p.d) * (t/((t-p.d)^2+1)))} ### add check statement in loop for function name
f.time.d.epi.w.var <- function(t,y){abs(y[1] + (y[2]*abs(y[3])/p.d) * (t/((t-p.d)^2+abs(y[3]))))}
f.time.d.cst <- function(t,y){rep(abs(y[1]),length(t))}
f.time.d.0 <- function(t,y){rep(0,length(t))}

vector.b.funcs.time <- c("b.cst"=f.time.b.cst,
                         "b.exp.abs"=f.time.b.exp.abs,
                         "b.exp.max"=f.time.b.exp.max,
                         "b.lin.abs"=f.time.b.lin.abs,
                         "b.lin.max"=f.time.b.lin.max,
                         "b.epi.fix"=f.time.b.epi.w.fix,
                         "b.epi.var"=f.time.b.epi.w.var)

vector.d.funcs.time <- c("d.0"=f.time.d.0,
                         "d.cst"=f.time.d.cst,
                         "d.exp.abs"=f.time.d.exp.abs,
                         "d.exp.max"=f.time.d.exp.max,
                         "d.lin.abs"=f.time.d.lin.abs,
                         "d.lin.max"=f.time.d.lin.max,
                         "d.epi.fix"=f.time.d.epi.w.fix,
                         "d.epi.var"=f.time.d.epi.w.var)

if(!("RPANDA_RTT_data" %in% list.files(".."))){
  system2("mkdir",args="../RPANDA_RTT_data")
}

for(i in 1:nrow(df)){
  # for new data
  func.name.df <- df$Best.model[i]
  if(func.name.df == "CB"){
    lamb.f.df <- vector.b.funcs.time[which(names(vector.b.funcs.time)=="b.cst")]
    mu.f.df <- vector.d.funcs.time[which(names(vector.d.funcs.time)=="d.0")]
  }else if(func.name.df == "CBD"){
    lamb.f.df <- vector.b.funcs.time[which(names(vector.b.funcs.time)=="b.cst")]
    mu.f.df <- vector.d.funcs.time[which(names(vector.d.funcs.time)=="d.cst")]
  }else{
    lamb.name.df <- substr(func.name.df,1,(unlist(gregexpr('_', func.name.df))[1]-1))
    mu.name.df <- substr(func.name.df,(unlist(gregexpr('_', func.name.df))[1]+1),nchar(func.name.df))
    
    if(substr(lamb.name.df,3,5) == "epi"){
      p.b <- as.numeric(substr(lamb.name.df,11,nchar(lamb.name.df)))
      lamb.name.df <- substr(lamb.name.df,1,9)
    }
    
    if(substr(mu.name.df,3,5) == "epi"){
      p.d <- as.numeric(substr(mu.name.df,11,nchar(mu.name.df)))
      mu.name.df <- substr(mu.name.df,1,9)
    }
    
    lamb.f.df <- vector.b.funcs.time[which(names(vector.b.funcs.time)==lamb.name.df)]
    mu.f.df <- vector.d.funcs.time[which(names(vector.d.funcs.time)==mu.name.df)]
  }
  
  timeseries <- timeseries_for_rates[,i]
  # you have to add omega and pi
  sp.rates.df <- as.data.frame(as.matrix(t(lamb.f.df[[1]](timeseries,as.numeric(c(df$lambda0[i],df$alpha[i],df$omega[i]))))))
  ex.rates.df <- as.data.frame(as.matrix(t(mu.f.df[[1]](timeseries,as.numeric(c(df$mu0[i],df$beta[i],df$pi[i]))))))
  net.div.rates.df <- sp.rates.df - ex.rates.df
  
  
  colnames(sp.rates.df) <- colnames(ex.rates.df) <- colnames(net.div.rates.df) <- timeseries
  
  write.csv(rbind(timeseries,sp.rates.df),paste0("../RPANDA_RTT_data/",df$Lineages[i],"_sp_rates_new.csv"))
  write.csv(rbind(timeseries,ex.rates.df),paste0("../RPANDA_RTT_data/",df$Lineages[i],"_ex_rates_new.csv"))
  write.csv(rbind(timeseries,net.div.rates.df),paste0("../RPANDA_RTT_data/",df$Lineages[i],"_nd_rates_new.csv"))
  
  
}

##########################################################################################################################

############################################# SCENARIO SELECTION #########################################################
criteria <- read.csv("../Scenarios_and_criteria.csv")
data <- read.csv("new_table_wo_env3.csv")
lin_scenarios <- matrix(ncol=2,nrow=nrow(data))
colnames(lin_scenarios) <- c('Lineage','Scenario')

ind.peak.dip <- grep("epi",data$best.model)
data.peak.dip <- data[ind.peak.dip,]
data <- data[-ind.peak.dip,]
data <- data[,-which(colnames(data) %in% c("omega","pi"))]

crit.colnames <- colnames(criteria)
crit <- unlist(criteria)
crit[which(crit=="0")] <- "==0"
criteria <- as.data.frame(matrix(crit,ncol = length(crit.colnames)))
colnames(criteria) <- crit.colnames

criteria.varying.but.gradual <- criteria[which(criteria$additional.remarks=="alpha=beta"),]
criteria.no.ex <- criteria[which(criteria$additional.remarks=="mu0=0"),]
criteria.3c <- criteria[which(criteria$additional.remarks=="alpha>beta"),]
criteria.no.extra <- criteria[which(criteria$additional.remarks==""),]


for(i in 1:nrow(data)){
  scenario <- c()
  best.params <- c(abs(data$lambda0[i]),data$alpha[i],abs(data$mu0[i]),data$beta[i])
  
  if(round(best.params[2],4)==round(best.params[4],4) & best.params[2]!=0 & best.params[4]!=0){
    criteria <- criteria.varying.but.gradual
  }else if(best.params[3] == 0){
    criteria <- criteria.no.ex
  }else if(best.params[2] > best.params[4]){
    criteria <- criteria.3c
  }else{
    criteria <- criteria.no.extra
  }
  params.comp <- c(best.params[2],best.params[4],(best.params[1]-best.params[3]))
  for(j in 1:nrow(criteria)){
    check <- parse(text=paste0(params.comp,criteria[j,4:6]))
    
    if(length(which(sapply(check,eval)))==3){
      scenario <- c(scenario,j)
    }
  }
  if(length(scenario)>1){
    scenario <- max(scenario)
  }
  
  if(is.null(scenario)){
    lin_scenarios[i,] <- c(data$lineage.names[i], "NA")
  }else{
    lin_scenarios[i,] <- c(data$lineage.names[i], criteria$Scenarios[scenario])
  }
  
}

ind.cri.na <- which(lin_scenarios[,2] == "NA")

# generally NA would be returned for cases where waxing and waning would be tricky
# hence just check from the RTT file and confirm if net-diversification went below 0 at some time point.
# then just assign SC 3 to it and leave the sub-scenario
for(i in ind.cri.na){
  nd.dat <- read.csv(paste0("../RPANDA_RTT_data/new_estimates/",lin_scenarios[i],"_nd_rates_new.csv"))
  if(length(which(nd.dat[2,] < 0))){
    lin_scenarios[i,2] <- "SC 3"
  }
  
}

for(i in 1:nrow(data.peak.dip)){
  t.f <- data.peak.dip[i,]$best.model
  
  t.f.b <- substr(t.f,1,(unlist(gregexpr('_', t.f))[1]-1))
  t.f.d <- substr(t.f,(unlist(gregexpr('_', t.f))[1]+1),nchar(t.f))
  
  pos.b.func <- pos.d.func <- NA
  
  if(nchar(t.f.b)>9){
    index.num <- unlist(lapply(0:9, gregexpr,t.f.b))
    index.num <- sort(index.num[-which(index.num == -1)])
    pos.b.func <-as.numeric(substr(t.f.b,index.num[1],nchar(t.f.b)))
    t.f.b <- substr(t.f.b,1,(index.num[1]-2))
  }
  
  if(nchar(t.f.d)>9){
    index.num <- unlist(lapply(0:9, gregexpr,t.f.d))
    index.num <- sort(index.num[-which(index.num == -1)])
    pos.d.func <- as.numeric(substr(t.f.d,index.num[1],nchar(t.f.d)))
    t.f.d <- substr(t.f.d,1,(index.num-2))
  }
  
  if(!is.na(pos.b.func) & !is.na(pos.d.func)){ # shift in both rates
    b.shift <- paste0(ifelse(data.peak.dip$alpha[i] > 0,"peak","dip"),"_b@",pos.b.func)
    d.shift <- paste0(ifelse(data.peak.dip$beta[i] > 0,"dip","peak"),"_d@",pos.d.func)
    lin_scenarios[nrow(data)+i,] <- c(data.peak.dip$lineage.names[i],paste0(b.shift,"+",d.shift))
  }else if(is.na(pos.b.func)){ # shift in just extinction rates
    if(data.peak.dip$alpha[i] == 0){ # + speciation rates are constant
      lin_scenarios[nrow(data)+i,] <- c(data.peak.dip$lineage.names[i],paste0(ifelse(data.peak.dip$beta[i] > 0,"dip","peak"),"_d@",pos.d.func))
    }else if(data.peak.dip$alpha[i] > 0){ # + speciation rates have decreased towards present
      lin_scenarios[nrow(data)+i,] <- c(data.peak.dip$lineage.names[i],paste0("SC 2 + ",ifelse(data.peak.dip$beta[i] > 0,"dip","peak"),"_d@",pos.d.func))
    }else{ # + speciation rates have increased towards present
      lin_scenarios[nrow(data)+i,] <- c(data.peak.dip$lineage.names[i],paste0("SC 4 + ",ifelse(data.peak.dip$beta[i] > 0,"dip","peak"),"_d@",pos.d.func))
    }
  }else{ # shift in just speciation rates
    if(data.peak.dip$beta[i] == 0){ # + extinction rates are constant
      lin_scenarios[nrow(data)+i,] <- c(data.peak.dip$lineage.names[i],paste0(ifelse(data.peak.dip$alpha[i] > 0,"peak","dip"),"_b@",pos.b.func))
    }else if(data.peak.dip$beta[i] > 0){ # + extinction rates have decreased towards present
      lin_scenarios[nrow(data)+i,] <- c(data.peak.dip$lineage.names[i],paste0("SC 4 + ",ifelse(data.peak.dip$alpha[i] > 0,"peak","dip"),"_b@",pos.d.func))
    }else{ # + speciation rates have increased towards present
      lin_scenarios[nrow(data)+i,] <- c(data.peak.dip$lineage.names[i],paste0("SC 2 + ",ifelse(data.peak.dip$alpha[i] > 0,"peak","dip"),"_b@",pos.d.func))
    }
  }
}

lin_scenarios <- lin_scenarios[order(lin_scenarios[,1],decreasing = FALSE),]

write.csv(lin_scenarios,"../model_comparisons_parallelized/lineages_scenarios.csv")
##############################################################################################################################

########################### SUMMARISE DRIVERS OF DIVERSIFICATION ##############################
#####################################################################################
setwd("../model_comparisons_parallelized/")
only.temp.mat <- read.csv("model_comparisons_w_smoothened_temp_model_concatenated_new.csv")
only.temp.mat <- cbind(only.temp.mat,matrix(nrow=nrow(only.temp.mat),ncol=12),"Lineage ID"=paste0(only.temp.mat$lineage.names,only.temp.mat$clades.node.nums))

him.temp.mat <- read.csv("model_comparisons_w_smoothened_him_temp_model_concatenated_very_new.csv")
him.temp.mat <- cbind(him.temp.mat,matrix(nrow=nrow(him.temp.mat),ncol=6),"Lineage ID"=paste0(him.temp.mat$lineage.names,him.temp.mat$clades.node.nums))

c4.him.temp.mat <- read.csv("model_comparisons_w_smoothened_c4_him_temp_model_concatenated_very_new.csv")
c4.him.temp.mat <- cbind(c4.him.temp.mat, "Lineage ID"=paste0(c4.him.temp.mat$lineage.names,c4.him.temp.mat$clades.node.nums))

colnames(only.temp.mat)
colnames(him.temp.mat)
colnames(c4.him.temp.mat)

# removing all the clades that are already there in the matrix that compares the effect of himalayan orogeny and temperature for avoiding repetition
only.temp.mat <- only.temp.mat[which(!only.temp.mat$`Lineage ID` %in% him.temp.mat$`Lineage ID`),]
# removing all the clades that are already there in the matrix that compares the effect of c4-expansion, himalayan orogeny and temperature for avoiding repetition
him.temp.mat <- him.temp.mat[which(!him.temp.mat$`Lineage ID` %in% c4.him.temp.mat$`Lineage ID`),]

# finally make a super matrix that has all the comparisons together
super.matrix <- as.data.frame(rbind(as.matrix(c4.him.temp.mat),as.matrix(him.temp.mat),as.matrix(only.temp.mat)))
super.matrix <- super.matrix[order(super.matrix$lineage.names,decreasing=FALSE),]
super.matrix <- super.matrix[,-1]
#View(super.matrix)

# code for selecting best model among the entire pool

# order it based on lineage names in alphabetical order
nms <- colnames(super.matrix) # names of columns of super.matrix

aicc.cols <- nms[c(which(grepl("AICc",nms)),which(grepl("aicc",nms)))] # extract just the columns that have AICc scores
aicc.cols <- aicc.cols[-which(aicc.cols == "Time.aicc")]

#aicc.super.matrix <- as.data.frame(cbind(super.matrix$lineage.names,super.matrix[,aicc.cols]))
super.matrix <- cbind(super.matrix,"lowest.aicc"=as.vector(apply(super.matrix[,aicc.cols], 1, min, na.rm = TRUE))) # add a column for the lowest AICc score per lineage
super.matrix <- cbind(super.matrix,"delta.aicc.CB"=(as.numeric(super.matrix[,aicc.cols[1]])-as.numeric(super.matrix[,"lowest.aicc"]))) # add a column for delta AICc for CB model
super.matrix <- cbind(super.matrix,"delta.aicc.CBD"=(as.numeric(super.matrix[,aicc.cols[2]])-as.numeric(super.matrix[,"lowest.aicc"]))) # add a column for delta AICc for CBD model
super.matrix <- cbind(super.matrix,"delta.aicc.DDD"=(as.numeric(super.matrix[,aicc.cols[3]])-as.numeric(super.matrix[,"lowest.aicc"]))) # # add a column for delta AICc for Time-varying model
super.matrix <- cbind(super.matrix,"delta.aicc.Temp"=(as.numeric(super.matrix[,aicc.cols[4]])-as.numeric(super.matrix[,"lowest.aicc"])))
super.matrix <- cbind(super.matrix,"delta.aicc.Him.oro"=(as.numeric(super.matrix[,aicc.cols[5]])-as.numeric(super.matrix[,"lowest.aicc"])))
super.matrix <- cbind(super.matrix,"delta.aicc.C4.Exp"=(as.numeric(super.matrix[,aicc.cols[6]])-as.numeric(super.matrix[,"lowest.aicc"])))

best.driver <- vector("character",nrow(super.matrix))

# Run through each lineage and check what's the best driver among - DDD, Temp, Him.oro and C4-Exp
for(i in 1:nrow(super.matrix)){
  drivers.w.eq.aicc <- as.vector(aicc.cols[which(super.matrix[i,(ncol(super.matrix)-5):ncol(super.matrix)]<2)])
  
  if(aicc.cols[1] %in% drivers.w.eq.aicc | aicc.cols[2] %in% drivers.w.eq.aicc){
    best.driver[i] <- 'Other'
    next
  }else {
    best.driver[i] <- paste(substr(drivers.w.eq.aicc,1,(nchar(drivers.w.eq.aicc)-5)),collapse = "+")
  }
  
}

super.matrix <- as.data.frame(cbind(super.matrix,"best.driver"=best.driver)) # add a column for the names of the best driver

#View(super.matrix)
write.csv(super.matrix,"../model_comparisons_parallelized/drivers_of_diversification.csv")

# collapse the conclusions of multiple sub-clades from same lineage and make a simpler matrix
# super.matrix <- read.csv("../model_comparisons_parallelized/drivers_of_diversification.csv")
lineages.w.drivers <- unique(super.matrix$lineage.names)
all.driver.combinations <- unique(super.matrix$best.driver[-which(super.matrix$best.driver == "Other")])
lineages <- read.csv("new_table_wo_env3.csv")$lineage.names
all.drivers <- matrix(nrow=length(lineages),ncol=6)
row.names(all)
colnames(all.drivers) <- c("Lineages","DDD","Temp","Him.oro","C4-Exp","Other")

for(i in 1:length(lineages)){
  #i=4
  ind.l <- which(super.matrix$lineage.names == lineages[i])
  drivers <- all.driver.combinations[all.driver.combinations  %in% unique(super.matrix$best.driver[ind.l])]
  if(length(drivers)==0){
    all.drivers[i,] <- c(lineages[i],0,0,0,0,1)
  }else{
    grep.DDD <- grep("DDD",drivers) %>% length %>% {ifelse(.>0,1,0)}
    grep.Temp <- grep("Temp",drivers) %>% length %>% {ifelse(.>0,1,0)}
    grep.Him.oro <- grep("Him.oro",drivers) %>% length %>% {ifelse(.>0,1,0)}
    grep.C4.exp <- grep("C4.Exp",drivers) %>% length %>% {ifelse(.>0,1,0)}
    all.drivers[i,] <- c(lineages[i],grep.DDD,grep.Temp,grep.Him.oro,grep.C4.exp,NA)
  }
}
lineages.not.present <- lineage.names[which(!lineage.names %in% lineages)]

all.drivers <- rbind(all.drivers,
                     cbind(lineages.not.present,matrix(nrow=length(lineages.not.present),ncol=4),
                           rep(1,length(lineages.not.present))))
all.drivers <- all.drivers[order(all.drivers[,1],decreasing=FALSE),]
write.csv(all.drivers,"../model_comparisons_parallelized/drivers-simplified.csv")

#######################################################################################################



#################################### MAKE FINAL SUMMARY ###################################################
###########################################################################################################
###################### PLEASE USE MS EXCEL/GOOGLE SHEETS/LIBREOFFICE CALC :-) #############################
# Just use results from "model_comparisons_wo_env.csv"(speciation rates, extinction rates and net-diversification rates under CBD)
# ,"lin_scenarios.csv" and "drivers-simplified.csv"

setwd("../../")
