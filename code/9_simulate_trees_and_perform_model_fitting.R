# simulations to validate RPANDA estimates

library(ape)
library(phytools)
library(RPANDA)
library(TESS)
library(ggplot2)
library(gridExtra)
library(grid)
library(Cairo)
library(dplyr)
library(tidyr)
library(ggbreak)

setwd("data/model_comparisons_parallelized/") # this folder would get created once code 3 is run
list.files()
dir.create("../Simulations")
dir.create("../Simulations/simulated_trees")
# read model parameter estimates from RPANDA
df <- read.csv("new_table_wo_env3.csv")

# TESS provides a ready-made function for simulating time-dependent BD
# but the interpretation of time in TESS is exactly the opposite of how RPANDA does
# RPANDA uses time from present in their functions whereas,
# TESS uses time since origin in their simulation function
# hence the function names and the function definitions have been reversed in the following section


# pool of functions used in the analyses
# Functions for speciation
f.time.b.lin.abs <- function(t,y.l){abs(y.l[1] + y.l[2] * t)} 
f.time.b.lin.max <- function(t,y.l){sapply((y.l[1] + y.l[2] * t),max,0)}
f.time.b.exp.abs <- function(t,y.l){abs(y.l[1] * exp(y.l[2] * t))}
f.time.b.exp.max <- function(t,y.l){sapply((y.l[1] * exp(y.l[2] * t)),max,0)}
f.time.b.epi.w.fix <- function(t,y.l){y.l[1] + (y.l[2]*1/p.b) * (t/((t-p.b)^2+1))}
f.time.b.epi.w.var <- function(t,y.l){y.l[1] + (y.l[2]*abs(y.l[3])/p.b) * (t/((t-p.b)^2+abs(y.l[3])))}
f.time.b.cst <- function(t,y.l){abs(y.l[1])}

# Fucntions for extinctions
f.time.d.lin.abs <- function(t,y.m){abs(y.m[1] + y.m[2] * t)}
f.time.d.lin.max <- function(t,y.m){sapply((y.m[1] + y.m[2] * t),max,0)}
f.time.d.exp.abs <- function(t,y.m){abs(y.m[1] * exp(y.m[2] * t))}
f.time.d.exp.max <- function(t,y.m){sapply((y.m[1] * exp(y.m[2] * t)),max,0)}
f.time.d.epi.w.fix <- function(t,y.m){y.m[1] + (y.m[2]*1/p.d) * (t/((t-p.d)^2+1))} ### add check statement in loop for function name
f.time.d.epi.w.var <- function(t,y.m){y.m[1] + (y.m[2]*abs(y.m[3])/p.d) * (t/((t-p.d)^2+abs(y.m[3])))}
f.time.d.cst <- function(t,y.m){abs(y.m[1])}
f.time.d.0 <- function(t,y.m){0}

vector.b.funcs.time <- c("b.cst"=f.time.b.cst,
                         "b.exp.abs"=f.time.b.exp.abs,
                         "b.exp.max"=f.time.b.exp.max,
                         "b.lin.abs"=f.time.b.lin.abs,
                         "b.lin.max"=f.time.b.lin.max,
                         "b.epi.fix"=f.time.b.epi.w.fix,
                         "b.epi.var"=f.time.b.epi.w.var)
#names(vector.b.funcs.time)[c(6,7)] <- paste0(names(vector.b.funcs.time)[c(6,7)],".",c(p.b.fix,p.b.var))

vector.d.funcs.time <- c("d.0"=f.time.d.0,
                         "d.cst"=f.time.d.cst,
                         "d.exp.abs"=f.time.d.exp.abs,
                         "d.exp.max"=f.time.d.exp.max,
                         "d.lin.abs"=f.time.d.lin.abs,
                         "d.lin.max"=f.time.d.lin.max,
                         "d.epi.fix"=f.time.d.epi.w.fix,
                         "d.epi.var"=f.time.d.epi.w.var)
#names(vector.d.funcs.time)[c(7,8)] <- paste0(names(vector.d.funcs.time)[c(7,8)],".",c(p.d.fix,p.d.var))


# finding best-fitting functions
# 1st vector for lambda functions
# 2nd vector for mu functions
best.b.funcs <- best.d.funcs <- names.b.funcs <- names.d.funcs <- t.f <- pos.b.funcs <- pos.d.funcs <- c()

# loop over all the lineages and get info about their best-fitting speciation and extinction functions
for(i in 1:nrow(df)){
  #i<-29
  cat("\nCurrently checking Lineage:",df$lineage.names[i],"\n")
  t.f <- c(t.f,df[i,"best.model"])
  if(substr(t.f[i],1,2)!='CB'){
    t.f.b <- substr(t.f[i],1,(unlist(gregexpr('_', t.f[i]))[1]-1))
    t.f.d <- substr(t.f[i],(unlist(gregexpr('_', t.f[i]))[1]+1),nchar(t.f[i]))
    
    if(nchar(t.f.b)>9){
      index.num <- unlist(lapply(0:9, gregexpr,t.f.b))
      index.num <- sort(index.num[-which(index.num == -1)])
      pos.b.funcs <- c(pos.b.funcs, as.numeric(substr(t.f.b,index.num[1],nchar(t.f.b))))
      t.f.b <- substr(t.f.b,1,(index.num[1]-2))
    }else pos.b.funcs <- c(pos.b.funcs, NA)
    
    if(nchar(t.f.d)>9){
      index.num <- unlist(lapply(0:9, gregexpr,t.f.d))
      index.num <- sort(index.num[-which(index.num == -1)])
      pos.d.funcs <- c(pos.d.funcs, as.numeric(substr(t.f.d,index.num[1],nchar(t.f.d))))
      t.f.d <- substr(t.f.d,1,(index.num-2))
    }else pos.d.funcs <- c(pos.d.funcs, NA)
    
    best.b.funcs <- c(best.b.funcs,vector.b.funcs.time[which(names(vector.b.funcs.time)==t.f.b)])
    best.d.funcs <- c(best.d.funcs,vector.d.funcs.time[which(names(vector.d.funcs.time)==t.f.d)])
    
  }else if(t.f[i] == 'CBD'){
    best.b.funcs <- c(best.b.funcs,vector.b.funcs.time[1])
    best.d.funcs <- c(best.d.funcs,vector.d.funcs.time[2])
    pos.b.funcs <- c(pos.b.funcs,NA)
    pos.d.funcs <- c(pos.d.funcs,NA)
  }else{
    best.b.funcs <- c(best.b.funcs,vector.b.funcs.time[1])
    best.d.funcs <- c(best.d.funcs,vector.d.funcs.time[1])
    pos.b.funcs <- c(pos.b.funcs,NA)
    pos.d.funcs <- c(pos.d.funcs,NA)
  }
  names.b.funcs <- c(names.b.funcs,df$lineage.names[i])
  names.d.funcs <- c(names.d.funcs,df$lineage.names[i])
  cat("best b.func is: ",names(best.b.funcs[i])," & ")
  cat("best d.func is",names(best.d.funcs[i]),"\n")
}
names(best.b.funcs) <- names.b.funcs
names(best.d.funcs) <- names.d.funcs

# check: names(best.b.funcs) == df$lineage.names
# check: names(best.d.funcs) == df$lineage.names



# Simulate 1000 trees with these functions and the best model parameters
# here y[1] refers to speciation rate at the crown age and not the speciation rate at present
# here y[1] refers to extinction rate at the crown age and not the extinction rate at present
coeff_var.gamma <- coeff_var.ntips <- mean.diff.gamma <- mean.diff.ntips <- sd.diff.gamma <- sd.diff.ntips <- p.gamma <- p.ntips <- c()
sim.gamma.matrix <- matrix(nrow=nrow(df),ncol=1000)
sim.ntips.matrix <- matrix(nrow=nrow(df),ncol=1000)

rownames(sim.gamma.matrix) <- rownames(sim.ntips.matrix) <- df$lineage.names

cr.ages <- df$Clade.age

for(i in 1:nrow(df)){
  #i=3
  lambda0.rpanda <- df$lambda0[i]
  mu0.rpanda <- df$mu0[i]
  al.rpanda <- df$alpha[i]
  be.rpanda <- df$beta[i]
  omega.rpanda <- df$omega[i]
  pi.rpanda <- df$pi[i]
  y.l <- c(lambda0.rpanda,al.rpanda,omega.rpanda)
  y.m <- c(mu0.rpanda,be.rpanda,pi.rpanda)
  p.b <- pos.b.funcs[i]
  p.d <- pos.d.funcs[i]
  
  lambda0.tess <- best.b.funcs[[i]](cr.ages[i],y.l)
  mu0.tess <- best.d.funcs[[i]](cr.ages[i],y.m)
  
  p.b <- cr.ages[i] - p.b
  p.d <- cr.ages[i] - p.d
  
  t.f.b <- substr(t.f[i],1,(unlist(gregexpr('_', t.f[i]))[1]-1))
  t.f.d <- substr(t.f[i],(unlist(gregexpr('_', t.f[i]))[1]+1),nchar(t.f[i]))
  
  if(substr(t.f.b,3,5)=="epi"){
    t.f.b <- substr(t.f.b,1,9)
    al.tess <- al.rpanda
  }else{
    al.tess <- -al.rpanda
  }
  
  if(substr(t.f.d,3,5)=="epi"){
    t.f.d <- substr(t.f.d,1,9)
    be.tess <- be.rpanda
  }else{
    be.tess <- -be.rpanda
  }
  omega.tess <- omega.rpanda
  pi.tess <- pi.rpanda
  
  y.l <- c(lambda0.tess,al.tess,omega.tess)
  y.m <- c(mu0.tess,be.tess,pi.tess)
  
  if(substr(t.f[i],1,2)=='CB'){
    b.func.tess <- function(t){lambda0.tess}
    d.func.tess <- function(t){mu0.tess}
  }else{
    # Functions for speciation
    f.time.b.lin.abs <- function(t){abs(y.l[1] + y.l[2] * t)} 
    f.time.b.lin.max <- function(t){sapply((y.l[1] + y.l[2] * t),max,0)}
    f.time.b.exp.abs <- function(t){abs(y.l[1] * exp(y.l[2] * t))}
    f.time.b.exp.max <- function(t){sapply((y.l[1] * exp(y.l[2] * t)),max,0)}
    f.time.b.epi.w.fix <- function(t){y.l[1] + (y.l[2]*1/p.b) * (t/((t-p.b)^2+1))}
    f.time.b.epi.w.var <- function(t){y.l[1] + (y.l[2]*abs(y.l[3])/p.b) * (t/((t-p.b)^2+abs(y.l[3])))}
    f.time.b.cst <- function(t){abs(y.l[1])}
    
    # Fucntions for extinctions
    f.time.d.lin.abs <- function(t){abs(y.m[1] + y.m[2] * t)}
    f.time.d.lin.max <- function(t){sapply((y.m[1] + y.m[2] * t),max,0)}
    f.time.d.exp.abs <- function(t){abs(y.m[1] * exp(y.m[2] * t))}
    f.time.d.exp.max <- function(t){sapply((y.m[1] * exp(y.m[2] * t)),max,0)}
    f.time.d.epi.w.fix <- function(t){y.m[1] + (y.m[2]*1/p.d) * (t/((t-p.d)^2+1))} ### add check statement in loop for function name
    f.time.d.epi.w.var <- function(t){y.m[1] + (y.m[2]*abs(y.m[3])/p.d) * (t/((t-p.d)^2+abs(y.m[3])))}
    f.time.d.cst <- function(t){abs(y.m[1])}
    f.time.d.0 <- function(t){0}
    
    vector.b.funcs.time <- c("b.cst"=f.time.b.cst,
                             "b.exp.abs"=f.time.b.exp.abs,
                             "b.exp.max"=f.time.b.exp.max,
                             "b.lin.abs"=f.time.b.lin.abs,
                             "b.lin.max"=f.time.b.lin.max,
                             "b.epi.fix"=f.time.b.epi.w.fix,
                             "b.epi.var"=f.time.b.epi.w.var)
    #names(vector.b.funcs.time)[c(6,7)] <- paste0(names(vector.b.funcs.time)[c(6,7)],".",c(p.b.fix,p.b.var))
    
    vector.d.funcs.time <- c("d.0"=f.time.d.0,
                             "d.cst"=f.time.d.cst,
                             "d.exp.abs"=f.time.d.exp.abs,
                             "d.exp.max"=f.time.d.exp.max,
                             "d.lin.abs"=f.time.d.lin.abs,
                             "d.lin.max"=f.time.d.lin.max,
                             "d.epi.fix"=f.time.d.epi.w.fix,
                             "d.epi.var"=f.time.d.epi.w.var)
    b.func.tess <- vector.b.funcs.time[[which(names(vector.b.funcs.time)==t.f.b)]] 
    d.func.tess <- vector.d.funcs.time[[which(names(vector.d.funcs.time)==t.f.d)]]
    
    
  }
  #trees <- tess.sim.age(n = 1000, age = cr.ages[i], lambda = b.func.tess, mu = d.func.tess)
  
  for(j in 1:100){
    #j=1
    #tr <- as.phylo(trees[[j]])
    chk=TRUE
    while(chk){
      tr <- as.phylo(tess.sim.age(n = 1, age = cr.ages[i], 
                                  lambda = b.func.tess, 
                                  mu = d.func.tess)[[1]])
      
      if(Ntip(tr)==df$N.tips[i]){chk=FALSE}
    }
    lineage_dir_path <- paste0("../Simulations/simulated_trees/",df$lineage.names[i])
    dir.create(lineage_dir_path)
    write.tree(phy = tr,file = paste0(lineage_dir_path,"/",df$lineage.names[i],"_",j,".tre"))
  }
}

############ RPANDA Model-fitting #####################
library(parallel)
library(foreach)
library(doParallel)
library(doFuture)

setwd("../Simulations/model_fits_of_simulated_trees/")
lapply(list.files("../simulated_trees/"),dir.create)
setwd("../simulated_trees/")
path_plus_filename <- function(name,pattern=NULL){
  return(paste0(name,"/",list.files(name,pattern=pattern)))
}

filenames <- unlist(lapply(list.files(),path_plus_filename))
# provide the tree file's name and the analysis would be performed
rpanda_models <- function(filename){
  # had to reload the packages as within the cluster made for parallelisation, everything outside of it becomes obsolete within it
  library(RPANDA)
  library(ape)
  library(phytools)
  library(phangorn)
  library(dispRity)
  library(dplyr)
  library(phangorn)
  library(TESS)
  library(Cairo)
  library(geoscale)
  library(unikn)
  
  # Estimate positions of potential peaks and dips in diversification rates (relevant to the function I designed to model symmetrical episodic rate shifts)
  estimate_peak_or_dip_position <- function(Tr,tot_time,res.const.bd,fix.width,peak.in.which.rate){
    res.epi <- c()
    
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
  # rpanda model
  # Function to fit and find the best parameter estimates for an time-varying and diversity-dependent diversification models  
  run_model_fit_wo_env <- function(Tr){
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
    
    
    if(Ntip(myTree)<4){
      cat("\nToo few tips to fit Time-varying models.\n")
      res.all <- as.data.frame(t(c("Clade.age"=tot_time,
                                   "CB.AICc"=res.const.b$aicc,
                                   "CB.lamb_par"=res.const.b$lamb_par,
                                   "CBD.AICc"=res.const.bd$aicc,
                                   "CBD.lamb_par"=res.const.bd$lamb_par,
                                   "CBD.mu_par"=res.const.bd$mu_par,
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
  
  model.comparisons.wo.env <- matrix(ncol = 16)
  colnames(model.comparisons.wo.env)<- c("lineage.names","N.tips","Clade.age","CB.AICc","CB.lamb_par","CBD.AICc","CBD.lamb_par","CBD.mu_par",
                                         "Best Time Model","Time.aicc","Time.lambda","Time.alpha","Time.omega","Time.mu","Time.beta","Time.pi")
  
  myTree <- read.tree(filename)
  
  
  cr.age <- max(branching.times(myTree))
  cat(paste0("\nCurrently using tree of ",substr(filename,1,(nchar(filename)-4)),"\n"))
  
  ################### NON-PALEOCLIMATIC_VARIABLE DEPENDENT MODELS ################
  # get the best paramters for CB, CBD, DDD and tD (time-dependence)
  res.npc1 <- run_model_fit_wo_env(myTree)
  
  # keep storing those values in this matrix 
  model.comparisons.wo.env <- rbind(model.comparisons.wo.env,
                                    c(substr(filename,1,(nchar(filename)-4)),
                                      Ntip(myTree),
                                      res.npc1$Clade.age,
                                      res.npc1$CB.AICc,
                                      res.npc1$CB.lamb_par,
                                      res.npc1$CBD.AICc,
                                      res.npc1$CBD.lamb_par,
                                      res.npc1$CBD.mu_par,
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
  write.csv(model.comparisons.wo.env[-1,],file = paste0("../model_fits_of_simulated_trees/",substr(filename,1,(nchar(filename)-4)),".csv"))
  
}

do <- function(pop, fun, ncores = 8, ...) {
  require(foreach)
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  foreach(i = pop) %dopar% fun(i, ...)
}
#

#
do(filenames, rpanda_models)

################## CoMET analysis ####################
system2("mkdir",args=c("../COMET_results_simulations"))
setwd("../COMET_results_simulations/")
lapply(list.files("../simulated_trees/"),dir.create)
system2("mkdir",args=c("../RTT_data_simulations"))
setwd("../RTT_data_simulations/")
lapply(list.files("../simulated_trees/"),dir.create)

setwd("../simulated_trees/")

perform_CoMET <- function(filename){
  
  library(TESS)
  myTree <- read.tree(filename)
  if(!is.ultrametric(myTree)) {myTree <- force.ultrametric(myTree,"extend")}
  cat("\nCurrently using tree of",substr(filename,1,(nchar(filename)-4)),"\n")
  times <- as.numeric(branching.times(myTree))
  cr.age <- max(times)
  ##### rjMCMC and CoMET analysis with empirical hyperpriors
  
  numExpectedMassExtinctions <- 1
  numExpectedRateChanges <- 1
  
  # lets keep it high
  expectedSurvivalProbability <- 0.5
  pMassExtinctionPriorShape2 <- 100
  pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 *
    expectedSurvivalProbability /
    (expectedSurvivalProbability - 1)
  
  samplingFraction <- 1
  
  # Make directory with same name as tree
  system2("mkdir",args=c(paste0("../COMET_results_simulations/",substr(filename,1,(nchar(filename)-4)))))
  setwd(paste0("../COMET_results_simulations/",substr(filename,1,(nchar(filename)-4))))
  
  # Perform TESS-COMET analysis
  convergence <- FALSE
  while(!convergence){
    tess.analysis(myTree,
                  empiricalHyperPriors = TRUE,
                  samplingProbability = samplingFraction,
                  estimateNumberMassExtinctions = FALSE,
                  MAX_ITERATIONS = 100000000,
                  dir = "COMET_without_mass_ex")
    
    
    output <- tess.process.output("COMET_without_mass_ex",
                                  numExpectedRateChanges = numExpectedRateChanges,
                                  numExpectedMassExtinctions = numExpectedMassExtinctions)
    sp.shift.df <- as.data.frame(output$`speciation shift times`)
    ex.shift.df <- as.data.frame(output$`extinction shift times`)
    
    bayes.factors <- as.data.frame(rbind(output$`speciation Bayes factors`,output$`extinction Bayes factors`))
    
    sp.rates.df <- as.data.frame(output$`speciation rates`)
    ex.rates.df <- as.data.frame(output$`extinction rates`)
    colnames(sp.rates.df) <- colnames(ex.rates.df) <- colnames(sp.shift.df) <- colnames(ex.shift.df) <- colnames(bayes.factors) <- output$intervals[-length(output$intervals)]
    
    
    
    #output$numExtinctionCategories
    #View(output$`net-diversification rates`)
    ess.vals <- c(
      effectiveSize(output$numSpeciationCategories),
      effectiveSize(output$numExtinctionCategories))
    
    if(length(which(ess.vals <= 200)) == 0){
      convergence <- TRUE
      
      write.csv(sp.rates.df,file=paste0("../../../RTT_data_simulations/",substr(filename,1,(nchar(filename)-4)),"_sp_rates.csv"))
      write.csv(ex.rates.df,file=paste0("../../../RTT_data_simulations/",substr(filename,1,(nchar(filename)-4)),"_ex_rates.csv"))
      write.csv(sp.shift.df,file=paste0("../../../RTT_data_simulations/",substr(filename,1,(nchar(filename)-4)),"_sp_shift_times.csv"))
      write.csv(ex.shift.df,file=paste0("../../../RTT_data_simulations/",substr(filename,1,(nchar(filename)-4)),"_ex_shift_times.csv"))
      write.csv(bayes.factors,file=paste0("../../../RTT_data_simulations/",substr(filename,1,(nchar(filename)-4)),"_bayes_factors.csv"))
    }
  }
  
  ### PLOTS TO CHECK FOR CONVERGENCE
  pdf("comet_no_mass_extinctions_convergence.pdf", 6, 6, bg="transparent")
  layout.mat <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
  layout(layout.mat)
  tess.plot.singlechain.diagnostics(output,
                                    parameters = c("speciation rates",
                                                   "extinction rates"),
                                    las=2)
  
  dev.off()
  
  
  
  ### PLOTS FOR RATES
  pdf("rates.pdf", 6, 6, bg="transparent")
  par(mfrow=c(2,2),omi = c(0,0,0.5,0.3))
  
  # Speciation rates
  tess.plot.output(output,
                   fig.types = c("speciation rates",
                                 "speciation shift times"),
                   las=2,
                   col = c("#8F00FF","#8F00FF"))
  
  # Extinction rates
  tess.plot.output(output,
                   fig.types = c("extinction rates",
                                 "extinction shift times"),
                   las=2,col=c("#FF0000","#FF0000"))
  
  
  dev.off()
  setwd("../../../simulated_trees/")
}

do(filenames,perform_CoMET)

########## Finding times of rate shifts - CoMET analysis ############
setwd("../RTT_data_simulations/")

# extract bayes factors through time bins
bayes_factor_files <- unlist(lapply(list.files(),path_plus_filename,"bayes"))

# find out times of strong rate-shifts
sp.rate.shift.times <- ex.rate.shift.times <- rep(NA, length(bayes_factor_files))
names(sp.rate.shift.times) <- names(ex.rate.shift.times) <- lineage.names

# loop through the time bins and find out the ones where bayes factor is >4.6 (i.e. support for episodic rate shift is high enough)
# check methods for justification of the threshold
for(i in 1:length(bayes_factor_files)){
  #i=23
  bf <- read.csv(bayes_factor_files[i])[,-1]
  sp.bf <- bf[1,]
  ex.bf <- bf[2,]
  times <- as.numeric(substr(colnames(bf),2,nchar(colnames(bf))))
  colnames(sp.bf) <- colnames(ex.bf) <- times
  
  high.sp.bf <- sp.bf[which(sp.bf>4.6)]
  high.ex.bf <- ex.bf[which(ex.bf>4.6)]
  sp.rate.shift.times[i] <- paste(colnames(high.sp.bf),collapse=",")
  ex.rate.shift.times[i] <- paste(colnames(high.ex.bf),collapse = ",")
  
}
# write a file that summaries the times of rate shifts per lineages 
write.csv(cbind("Lineage names"=bayes_factor_files,
                "Speciation rate shift times"=as.vector(sp.rate.shift.times),
                "Extinction rate shift times"=as.vector(ex.rate.shift.times)),"../times_of_strong_rate_shifts.csv")

############# Scenario-assessment, stats for estimations ##############


assess_scenario <- function(model_fitting_matrix){
  # code for selecting best model among the entire pool
  df <- model_fitting_matrix
  #View(df)
  
  # order it based on lineage names in alphabetical order
  df <- df[order(df$lineage.names,decreasing=FALSE),]
  nms <- colnames(df) # names of columns of df
  
  
  
  aicc.cols <- nms[c(which(grepl("AICc",nms)),which(grepl("aicc",nms)))] # extract just the columns that have AICc scores
  
  df <- cbind(df,"lowest.aicc"=as.vector(apply(df[,aicc.cols], 1, min, na.rm = TRUE))) # add a column for the lowest AICc score per lineage
  df <- cbind(df,"delta.aicc.CB"=(as.numeric(df[,aicc.cols[1]]) - as.numeric(df[,"lowest.aicc"]))) # add a column for delta AICc for CB model
  df <- cbind(df,"delta.aicc.CBD"=(as.numeric(df[,aicc.cols[2]]) - as.numeric(df[,"lowest.aicc"]))) # add a column for delta AICc for CBD model
  df <- cbind(df,"delta.aicc.Time"=(as.numeric(df[,aicc.cols[3]])-as.numeric(df[,"lowest.aicc"]))) # # add a column for delta AICc for Time-varying model
  
  best.model <- vector("character",nrow(df))
  best.model.paramters <- matrix(nrow=1,ncol=6)
  colnames(best.model.paramters) <- c("lambda0","alpha","omega","mu0","beta","pi")
  
  # Run through each lineage and check what's the best model among - CB, CBD and Time-varying
  for(i in 1:nrow(df)){
    #i=1
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
  if(best.model.paramters["omega"] == "0") best.model.paramters["omega"] <- "NA"
  if(best.model.paramters["pi"] == "0") best.model.paramters["pi"] <- "NA"
  
  cat("\n Best model parameters are: ",best.model.paramters,"\n")
  
  df <- as.data.frame(cbind(df,"best.model"=best.model,as.data.frame(t(best.model.paramters)))) # add a column for the names of the best model out of : CB, CBD, Time
  ind.time.best <- which(df$best.model=="Time") # get the indices for all the lineages/rows where "Time" is the best model
  df[ind.time.best,"best.model"] <- df$`Best Time Model`[ind.time.best] # replace "Time" with the actual time-varying model name in this section of the matrix
  
  # SCENARIO ASSESSMENT
  criteria <- read.csv("../../Scenarios_and_criteria.csv")
  data <- df
  ind.peak.dip <- grep("epi",data$best.model)
  
  lin_scenarios <- matrix(ncol=2,nrow=nrow(data))
  colnames(lin_scenarios) <- c('Lineage','Scenario')
  
  if(length(ind.peak.dip)>0){
    # data <- data[-ind.peak.dip,]
    # data <- data[,-which(colnames(data) %in% c("omega","pi"))]
    
    t.f <- data$best.model
    
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
      b.shift <- paste0(ifelse(data$alpha > 0,"peak","dip"),"_b@",pos.b.func)
      d.shift <- paste0(ifelse(data$beta > 0,"dip","peak"),"_d@",pos.d.func)
      lin_scenarios <- c(data$lineage.names,paste0(b.shift,"+",d.shift))
    }else if(is.na(pos.b.func)){ # shift in just extinction rates
      if(data$alpha == 0){ # + speciation rates are constant
        lin_scenarios <- c(data$lineage.names,paste0(ifelse(data$beta > 0,"dip","peak"),"_d@",pos.d.func))
      }else if(data$alpha > 0){ # + speciation rates have decreased towards present
        lin_scenarios <- c(data$lineage.names,paste0("SC 2 + ",ifelse(data$beta > 0,"dip","peak"),"_d@",pos.d.func))
      }else{ # + speciation rates have increased towards present
        lin_scenarios <- c(data$lineage.names,paste0("SC 4 + ",ifelse(data$beta > 0,"dip","peak"),"_d@",pos.d.func))
      }
    }else{ # shift in just speciation rates
      if(data$beta == 0){ # + extinction rates are constant
        lin_scenarios <- c(data$lineage.names,paste0(ifelse(data$alpha > 0,"peak","dip"),"_b@",pos.b.func))
      }else if(data$beta > 0){ # + extinction rates have decreased towards present
        lin_scenarios <- c(data$lineage.names,paste0("SC 4 + ",ifelse(data$alpha > 0,"peak","dip"),"_b@",pos.b.func))
      }else{ # + speciation rates have increased towards present
        lin_scenarios <- c(data$lineage.names,paste0("SC 2 + ",ifelse(data$alpha > 0,"peak","dip"),"_b@",pos.b.func))
      }
    }
    
  }else{
    crit.colnames <- colnames(criteria)
    crit <- unlist(criteria)
    crit[which(crit=="0")] <- "==0"
    criteria <- as.data.frame(matrix(crit,ncol = length(crit.colnames)))
    colnames(criteria) <- crit.colnames
    
    criteria.varying.but.gradual <- criteria[which(criteria$additional.remarks=="alpha=beta"),]
    criteria.no.ex <- criteria[which(criteria$additional.remarks=="mu0=0"),]
    criteria.3c <- criteria[which(criteria$additional.remarks=="alpha>beta"),]
    criteria.no.extra <- criteria[which(criteria$additional.remarks==""),]
    
    
    #i=16
    #data$lineage.names[16]
    scenario <- c()
    best.params <- c(abs(as.numeric(data$lambda0[i])),
                     as.numeric(data$alpha[i]),
                     abs(as.numeric(data$mu0[i])),
                     as.numeric(data$beta[i]))
    
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
      #j=1
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
    
    # generally NA would be returned for cases where waxing and waning would be tricky
    # hence just check from the RTT file and confirm if net-diversification went below 0 at some time point.
    # then just assign SC 3 to it and leave the sub-scenario
    if(length(which(lin_scenarios[,2] == "NA")) != 0){
      cat("it turned out to be NA. lineage:",data$lineage.names,"\n")
      sp.dat <- read.csv(paste0("../RTT_data_simulations/",data$lineage.names,"_sp_rates.csv"))
      ex.dat <- read.csv(paste0("../RTT_data_simulations/",data$lineage.names,"_ex_rates.csv"))
      nd.dat <- sp.dat[,-1] - ex.dat[,-1]
      if(length(which(nd.dat < 0))>0){
        lin_scenarios[1,2] <- "SC 3"
      }
    }
    
  }
  return(list("scenarios"=lin_scenarios,"Best.model.parameters"=best.model.paramters))
}

compute_stats_for_parameters <- function(empirical_val,estimates_distribution){
  # w1.sum <- wilcox.test(estimates_distributio, mu = empirical_val, alternative = "two.sided")
  # p <- w1.sum$p.value
  
  return(estimates_distribution - empirical_val)/sd(estimates_distribution)
  # mean.diff <- c(mean.diff,mean(norm.diff))
  # sd.diff <- c(sd.diff,sd(norm.diff))
  #coeff_var <- c(coeff_var,(sd.diff[i]/mean[i]))
}

setwd("../model_fits_of_simulated_trees/")
folders <- list.files()

empirical_scenarios <- read.csv("../../Codes_and_Data/data/model_comparisons_parallelized/lineages_scenarios.csv")
empirical_best_model_parameters <- read.csv("../../Codes_and_Data/data/model_comparisons_parallelized/new_table_wo_env3.csv")

scenario_estimation_proportions <- matrix(nrow=length(folders),ncol=7)
sim.lambda0.matrix <- matrix(nrow=length(folders),ncol=100)
sim.mu0.matrix <- matrix(nrow=length(folders),ncol=100)
sim.alpha.matrix <- matrix(nrow=length(folders),ncol=100)
sim.beta.matrix <- matrix(nrow=length(folders),ncol=100)
sim.b.shift.matrix <- matrix(nrow=length(folders),ncol=100)
sim.d.shift.matrix <- matrix(nrow=length(folders),ncol=100)

rownames(sim.lambda0.matrix) <- rownames(sim.mu0.matrix) <- rownames(sim.alpha.matrix) <- 
  rownames(sim.beta.matrix) <- rownames(sim.b.shift.matrix) <- rownames(sim.d.shift.matrix) <- folders
colnames(scenario_estimation_proportions) <- c("Lineage","Ntips","Clade_age","#NAs","Prop_Correct_Inference","Prop_Incorrect_Inference","Empirical Scenario")
# loop through all the model parameters and assess scenarios
for(i in 1:length(folders)){
  #i <- 3
  scenario_for_this_lineage <- refined_scenario_for_this_lineage <- c()
  
  csv_files <- list.files(folders[i])
  num_files <- length(csv_files)
  for(j in 1:num_files){
    #j <- 14
    model_fit_data <- as.data.frame(t(read.csv(paste0("./",folders[i],"/",csv_files[j]))))
    colnames(model_fit_data) <- model_fit_data[1,]
    model_fit_data <- model_fit_data[-1,]
    scenario_and_best_model_parameters <- assess_scenario(model_fitting_matrix = model_fit_data)
    just_scenario <- scenario_and_best_model_parameters$scenarios[2]
    
    
    sim.lambda0.matrix[i,j] <- abs(as.numeric(scenario_and_best_model_parameters$Best.model.parameters["lambda0"]))
    sim.mu0.matrix[i,j] <- abs(as.numeric(scenario_and_best_model_parameters$Best.model.parameters["mu0"]))
    sim.alpha.matrix[i,j] <- as.numeric(scenario_and_best_model_parameters$Best.model.parameters["alpha"])
    sim.beta.matrix[i,j] <- as.numeric(scenario_and_best_model_parameters$Best.model.parameters["beta"])
    
    ind.at <- unlist(gregexpr("@",just_scenario))
    
    sim.b.shift.matrix[i,j] <- ifelse(length(grep("b@",just_scenario))!=0,
                                      as.numeric(substr(just_scenario,
                                                        unlist(gregexpr("b@",just_scenario))[1]+2,
                                                        ifelse(length(ind.at) > 1,(unlist(gregexpr("\\+",just_scenario))[1]-1),nchar(just_scenario))
                                      )),
                                      0)
    sim.d.shift.matrix[i,j] <- ifelse(length(grep("d@",just_scenario))!=0,
                                      as.numeric(substr(just_scenario,
                                                        unlist(gregexpr("d@",just_scenario))[1]+2,
                                                        nchar(just_scenario))),
                                      0)
    if(length(ind.at)==1){
      if(ind.at!=-1){
        just_scenario <- substr(just_scenario,1,ind.at-1)
      }
    }else if(length(ind.at) > 1){
      part_b <- substr(just_scenario,1,ind.at[1]-1)
      part_d <- substr(just_scenario,unlist(gregexpr("\\+",just_scenario))[1]+1,ind.at[2]-1)
      just_scenario <- paste0(part_b,"+",part_d)
    }
    scenario_for_this_lineage <- c(scenario_for_this_lineage,
                                   just_scenario)
  }
  emp.scenario <- empirical_scenarios$Scenario[i]
  ind.at <- unlist(gregexpr("@",emp.scenario))
  if(length(ind.at)==1){
    if(ind.at!=-1){
      emp.scenario <- substr(emp.scenario,1,ind.at-1)
    }
  }else if(length(ind.at) > 1){
    part_b <- substr(emp.scenario,1,ind.at[1]-1)
    part_d <- substr(emp.scenario,unlist(gregexpr("\\+",emp.scenario))[1]+1,ind.at[2]-1)
    emp.scenario <- paste0(part_b,"+",part_d)
  }
  
  # if(empirical_scenarios$Scenario[i] == "SC 3"){
  #   prop_correct <- length(which(substr(scenario_for_this_lineage,1,4) == "SC 3"))/num_files
  # }
  
  
  prop_correct <- length(which(unlist(gregexpr(substr(emp.scenario,1,4),scenario_for_this_lineage))!=-1))/num_files
  scenario_estimation_proportions[i,] <- c(folders[i],
                                           read.csv(paste0("./",folders[i],"/",csv_files[1]))$x[c(2,3)],
                                           length(which(scenario_for_this_lineage=="NA")),
                                           prop_correct,
                                           (1-prop_correct),
                                           empirical_scenarios$Scenario[i])
  sim.lambda0.matrix[i,] <- compute_stats_for_parameters(empirical_best_model_parameters$lambda0[i],sim.lambda0.matrix[i,])
  sim.mu0.matrix[i,] <- compute_stats_for_parameters(empirical_best_model_parameters$mu0[i],sim.mu0.matrix[i,])
  sim.alpha.matrix[i,] <- compute_stats_for_parameters(empirical_best_model_parameters$alpha[i],sim.alpha.matrix[i,])
  sim.beta.matrix[i,] <- compute_stats_for_parameters(empirical_best_model_parameters$beta[i],sim.beta.matrix[i,])
  
  emp.best.model.name <- empirical_best_model_parameters$best.model[i]
  b.shift <- d.shift <- 0
  if(length(grep("b.epi",emp.best.model.name))>0){
    b.shift <- as.numeric(substr(emp.best.model.name,
                                 unlist(gregexpr("\\.",emp.best.model.name))[3]+1,
                                 unlist(gregexpr("_",emp.best.model.name))-1))
  }else if(length(grep("d.epi",emp.best.model.name))>0){
    d_part <- substr(emp.best.model.name,unlist(gregexpr("_",emp.best.model.name))+1,nchar(emp.best.model.name))
    d.shift <- as.numeric(substr(d_part,
                                 unlist(gregexpr("\\.",emp.best.model.name))[3]+1,
                                 unlist(gregexpr("_",emp.best.model.name))-1))
  }
  #sim.b.shift.matrix[i,] <- compute_stats_for_parameters(b.shift,sim.b.shift.matrix[i,])
  #sim.d.shift.matrix[i,] <- compute_stats_for_parameters(d.shift,sim.d.shift.matrix[i,])
}

scenario_estimation_proportions <- as.data.frame(scenario_estimation_proportions)
scenario_estimation_proportions$Ntips <- as.numeric(scenario_estimation_proportions$Ntips)
scenario_estimation_proportions$Clade_age <- as.numeric(scenario_estimation_proportions$Clade_age)
scenario_estimation_proportions$Prop_Correct_Inference <- as.numeric(scenario_estimation_proportions$Prop_Correct_Inference)
scenario_estimation_proportions$Prop_Incorrect_Inference <- as.numeric(scenario_estimation_proportions$Prop_Incorrect_Inference)
scenario_estimation_proportions$`#NAs` <- as.numeric(scenario_estimation_proportions$`#NAs`)
scenario_estimation_proportions$Lineage[which(scenario_estimation_proportions$Prop_Correct_Inference < 0.5)]

scenario_estimation_proportions <- scenario_estimation_proportions[order(scenario_estimation_proportions$Ntips,decreasing = F),]
scenario_estimation_proportions$Lineage <- factor(scenario_estimation_proportions$Lineage,levels = scenario_estimation_proportions$Lineage)
View(scenario_estimation_proportions)


############### Visualisation ################
data <- read.csv("../../Codes_and_Data/data/Final_summary_updated.csv")
proportions_scenarios_long <- scenario_estimation_proportions %>% left_join(.,data %>% dplyr::select(Taxonomic.group, Lineage) %>% distinct()) %>% 
  arrange(.,Ntips) %>% mutate(.,name_clade_size=paste0(Lineage,"(",Ntips,")")) %>% 
  dplyr::select(.,name_clade_size,Taxonomic.group,Prop_Incorrect_Inference,Prop_Correct_Inference) %>%
  pivot_longer(cols=starts_with("Pr"),names_to="category",values_to = "Proportions")

proportions_scenarios_long$category <- factor(proportions_scenarios_long$category,levels = c("Prop_Incorrect_Inference","Prop_Correct_Inference"))

prop_correct_vs_lineages <- ggplot(proportions_scenarios_long,aes(x=name_clade_size,y=Proportions,fill=category))+
  geom_bar(lwd=0.2,stat="identity",color="black")+
  scale_fill_manual(values=c("white","black"),labels=c("Incorrect","Correct"))+
  facet_grid(facets = vars(Taxonomic.group)) +
  ylab("Proportion of correctly inferred diversification scenarios")+
  xlab("Lineages")+
  theme_bw()+
  ggtitle("Diversification Scenario prediction across lineages of varying tree size")+
  theme(legend.position = "top", axis.text.x=element_text(angle=45,hjust=1),plot.title = element_text(hjust = 0.5))

scenario_estimation_proportions$simplified_scenario <- scenario_estimation_proportions$`Empirical Scenario`
scenario_estimation_proportions$simplified_scenario[grep("@",scenario_estimation_proportions$simplified_scenario)] <- "Shift"
scenario_estimation_proportions$simplified_scenario[grep("SC",scenario_estimation_proportions$simplified_scenario)] <- substr(scenario_estimation_proportions$simplified_scenario[grep("SC",scenario_estimation_proportions$simplified_scenario)],1,4)

prop_correct_vs_scenarios <- ggplot(scenario_estimation_proportions,aes(x=simplified_scenario,y=Prop_Correct_Inference))+
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2),size=3)+
  xlab("Diversification scenarios")+ylab("Proportion of correctly inferred diversification scenarios")+theme_bw()

prop_correct_vs_ntips <- ggplot(scenario_estimation_proportions,aes(x=Ntips,y=Prop_Correct_Inference))+
  geom_point(shape=16, size=3)+
  xlab("No. of tips")+ylab("Proportion of correctly inferred diversification scenarios")+theme_bw()

pdf("../scenario_estimations_RPANDA_models.pdf",11.7,8.3,bg="transparent")
grid.arrange(prop_correct_vs_lineages,arrangeGrob(prop_correct_vs_ntips,prop_correct_vs_scenarios,nrow=2),ncol=2)
dev.off()




lambda0 <- as.data.frame(t(sim.lambda0.matrix)) %>% pivot_longer(everything()) %>% dplyr::rename(Lineage=name)  %>% 
  left_join(.,data %>% dplyr::select(Taxonomic.group, Lineage) %>% distinct()) %>%
  ggplot(aes(x=Lineage,y=value))+#facet_grid(facets = vars(Taxonomic.group)) +
  ylab("Standard effect size(\u03bb0)")+
  geom_hline(yintercept = c(-2,2), linetype = 2,size=0.2,color="red")+
  geom_boxplot(lwd=0.1,outlier.size = 0.5)+theme_bw()+
  #ggtitle("Standard effect size: comparing empirical and 100 simulated trees")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(angle = 90))

mu0 <- as.data.frame(t(sim.mu0.matrix)) %>% pivot_longer(everything()) %>% dplyr::rename(Lineage=name) %>% 
  left_join(.,data %>% dplyr::select(Taxonomic.group, Lineage) %>% distinct()) %>%
  ggplot(aes(x=Lineage,y=value))+#facet_grid(facets = vars(Taxonomic.group)) +
  ylab("Standard effect size(\u03bc0)")+
  geom_hline(yintercept = c(-2,2), linetype = 2,size=0.2,color="red")+
  geom_boxplot(lwd=0.1,outlier.size = 0.5)+theme_bw()+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(angle = 90))

alpha <- as.data.frame(t(sim.alpha.matrix)) %>% pivot_longer(everything()) %>% dplyr::rename(Lineage=name) %>% 
  left_join(.,data %>% dplyr::select(Taxonomic.group, Lineage) %>% distinct()) %>%
  ggplot(aes(x=Lineage,y=value))+#facet_grid(facets = vars(Taxonomic.group)) +
  ylab("Standard effect size(\u03b1)")+
  #scale_y_continuous(limits = c(-650,3))+
  scale_y_break(c(-3, -620),scale=2)+
  geom_hline(yintercept = c(-2,2), linetype = 2,size=0.2,color="red")+
  geom_boxplot(lwd=0.1,outlier.size = 0.5)+theme_bw()+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(angle = 90))

beta <- as.data.frame(t(sim.beta.matrix)) %>% pivot_longer(everything()) %>% dplyr::rename(Lineage=name) %>% 
  left_join(.,data %>% dplyr::select(Taxonomic.group, Lineage) %>% distinct()) %>%
  ggplot(aes(x=Lineage,y=value))+#facet_grid(facets = vars(Taxonomic.group)) +
  ylab("Standard effect size(\u03b2)")+
  #scale_y_continuous(limits = c(-3,3))+
  geom_hline(yintercept = c(-2,2), linetype = 2,size=0.2,color="red")+
  geom_boxplot(lwd=0.1,outlier.size = 0.5)+theme_bw()+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(angle = 90))
#theme(axis.text.x=element_text(angle=45,hjust=1),axis.text.y = element_text(angle = 90))
names.plot <-  as.data.frame(t(sim.beta.matrix)) %>% pivot_longer(everything()) %>% dplyr::rename(Lineage=name) %>% 
  left_join(.,data %>% dplyr::select(Taxonomic.group, Lineage) %>% distinct()) %>%
  ggplot(aes(x=Lineage,y=value))+#facet_grid(facets = vars(Taxonomic.group)) +
  ylab("Standard effect size(\u03b2)")+
  #scale_y_continuous(limits = c(-3,3))+
  #geom_hline(yintercept = c(-2,2), linetype = 2,size=0.2,color="red")+
  #geom_boxplot(lwd=0.1,outlier.size = 0.5)+
  theme_minimal()+
  theme(
    axis.title = element_blank(),       # Remove axis titles
    axis.text.y = element_blank(),       # Remove y-axis text
    axis.ticks.y = element_blank(),      # Remove y-axis ticks
    axis.line.y = element_blank(),       # Remove y-axis line
    #axis.ticks.x = element_blank(),      # Remove x-axis ticks
    axis.line.x = element_blank(),       # Remove x-axis line
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x=element_text(angle=45,hjust=1),
    plot.margin =  unit(c(1,0.8,1,0.8), "cm")# Remove gridlines
  )

Cairo::CairoPDF("../Model_parameter_estimations.pdf",8.3,11.7,bg="transparent")
grid.arrange(lambda0,mu0,alpha,beta,names.plot,nrow=5)
dev.off()


# comet_shift_times_data$Speciation.rate.shift.times
comet.b.shift.times <- comet.d.shift.times <- matrix(nrow=1,ncol=2)
colnames(comet.b.shift.times) <- colnames(comet.d.shift.times) <- c("Lineage","Time")

comet_shift_times_data <- read.csv("../times_of_strong_rate_shifts.csv")
head(comet_shift_times_data)

scenario_estimation_proportions$COMET.presence.of.shifts <- rep(NA,length(folders))
for(i in 1:length(folders)){
  #i <-1
  indices.rows <- grep(folders[i],comet_shift_times_data$Lineage.names)
  speciation.rate.shift.times.in.this.lineage <- comet_shift_times_data$Speciation.rate.shift.times[indices.rows]
  ind.trees.with.shifts.in.this.lineage <- which(speciation.rate.shift.times.in.this.lineage != "")
  scenario_estimation_proportions$COMET.presence.of.shifts[which(scenario_estimation_proportions$Lineage==folders[i])] <- length(ind.trees.with.shifts.in.this.lineage)
  
  for(j in ind.trees.with.shifts.in.this.lineage){
    #j <- ind.trees.with.shifts.in.this.lineage[2]
    if(length(grep(",",speciation.rate.shift.times.in.this.lineage[j]))==0){
      comet.b.shift.times <- rbind(comet.b.shift.times,c(folders[i], as.numeric(speciation.rate.shift.times.in.this.lineage[j])))
    }else{
      indices.of.comma <- c(0, # adding 1 here for convenience in cutting string
                            unlist(gregexpr(",",speciation.rate.shift.times.in.this.lineage[j])),
                            nchar(speciation.rate.shift.times.in.this.lineage[j])+1) # adding the index of the last character here for convenience in cutting string
      
      for(k in 1:(length(indices.of.comma)-1)){
        comet.b.shift.times <- rbind(comet.b.shift.times,
                                     c(folders[i], as.numeric(substr(speciation.rate.shift.times.in.this.lineage[j],indices.of.comma[k]+1,indices.of.comma[k+1]-1))))
      }
      
    }
    
  }
  
}
comet.b.shift.times <- comet.b.shift.times[-1,]


ind.of.lineages.with.emp.rate.shifts <- match(which(!(data$Speciation.rate.shift.times..TESS. %in% "")),as.numeric(rownames(scenario_estimation_proportions)))
scenario_estimation_proportions$prop.correct.COMET.presence.of.shifts <- rep(NA,length(folders))
scenario_estimation_proportions$prop.correct.COMET.presence.of.shifts[ind.of.lineages.with.emp.rate.shifts] <- (scenario_estimation_proportions$COMET.presence.of.shifts[ind.of.lineages.with.emp.rate.shifts])/100
scenario_estimation_proportions$prop.correct.COMET.presence.of.shifts[which(!(1:33 %in% ind.of.lineages.with.emp.rate.shifts))] <- (100 - scenario_estimation_proportions$COMET.presence.of.shifts[which(!(1:33 %in% ind.of.lineages.with.emp.rate.shifts))])/100
scenario_estimation_proportions$prop.incorrect.COMET.presence.of.shifts <- 1 - scenario_estimation_proportions$prop.correct.COMET.presence.of.shifts

proportions_scenarios_long <- scenario_estimation_proportions %>% left_join(.,data %>% dplyr::select(Taxonomic.group, Lineage) %>% distinct()) %>% 
  arrange(.,Ntips) %>% mutate(.,name_clade_size=paste0(Lineage,"(",Ntips,")")) %>% 
  dplyr::select(.,name_clade_size,Taxonomic.group,prop.correct.COMET.presence.of.shifts,prop.incorrect.COMET.presence.of.shifts) %>%
  pivot_longer(cols=starts_with("pr"),names_to="category",values_to = "Proportions")

proportions_scenarios_long$category <- factor(proportions_scenarios_long$category,levels = c("prop.incorrect.COMET.presence.of.shifts","prop.correct.COMET.presence.of.shifts"))

prop_correct_vs_lineages <- ggplot(proportions_scenarios_long,aes(x=name_clade_size,y=Proportions,fill=category))+
  geom_bar(lwd=0.2,stat="identity",color="black")+
  scale_fill_manual(values=c("white","black"),labels=c("Incorrect","Correct"))+
  facet_grid(facets = vars(Taxonomic.group)) +
  ylab("Proportion of correctly inferred diversification scenarios")+
  xlab("Lineages")+
  theme_bw()+
  ggtitle("Diversification Scenario prediction across lineages of varying tree size")+
  theme(legend.position = "top", axis.text.x=element_text(angle=45,hjust=1),plot.title = element_text(hjust = 0.5))

scenario_estimation_proportions$simplified_scenario <- scenario_estimation_proportions$`Empirical Scenario`
scenario_estimation_proportions$simplified_scenario[grep("@",scenario_estimation_proportions$simplified_scenario)] <- "Shift"
scenario_estimation_proportions$simplified_scenario[grep("SC",scenario_estimation_proportions$simplified_scenario)] <- substr(scenario_estimation_proportions$simplified_scenario[grep("SC",scenario_estimation_proportions$simplified_scenario)],1,4)

prop_correct_vs_scenarios <- ggplot(scenario_estimation_proportions,aes(x=simplified_scenario,y=prop.correct.COMET.presence.of.shifts))+
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2),size=3)+
  xlab("Diversification scenarios")+ylab("Proportion of correctly inferred diversification scenarios")+theme_bw()

prop_correct_vs_ntips <- ggplot(scenario_estimation_proportions,aes(x=Ntips,y=prop.correct.COMET.presence.of.shifts))+
  geom_point(shape=16, size=3)+
  xlab("No. of tips")+ylab("Proportion of correctly inferred diversification scenarios")+theme_bw()

pdf("../scenario_estimationsCoMET.pdf",11.7,8.3,bg="transparent")
grid.arrange(prop_correct_vs_lineages,arrangeGrob(prop_correct_vs_ntips,prop_correct_vs_scenarios,nrow=2),ncol=2)
dev.off()


comet.b.shift.times <- as.data.frame(comet.b.shift.times)
comet.b.shift.times$Time <- as.numeric(comet.b.shift.times$Time)

sim.b.shift.matrix <- as.data.frame(sim.b.shift.matrix)
sim.d.shift.matrix <- as.data.frame(sim.d.shift.matrix)
sim.b.shift.matrix[,101] <- rownames(sim.b.shift.matrix)
sim.d.shift.matrix[,101] <- rownames(sim.d.shift.matrix)

pdf("../times_of_strong_shift.pdf",11.7,8.3,bg="transparent")
b_shift_rpanda <- sim.b.shift.matrix %>% pivot_longer(cols = -V101,values_to = "Time") %>% dplyr::rename(Lineage=V101) %>% 
  ggplot(aes(x=Lineage,y=Time))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  ylim(c(0.5,160))+
  theme_bw()+
  geom_point(data=data, aes(y=Up_Peak_or_Down_Dip_Time),fill="blue",shape=24,color="white",size=3.5,alpha=0.7)+
  geom_point(data=data, aes(y=Peak_or_Dip_Time),fill="blue",shape=21,color="white",size=3.5,alpha=0.7)+
  geom_point(data=data, aes(y=Down_Peak_or_Up_Dip_Time),fill="blue",shape=25,color="white",size=3.5,alpha=0.7)+
  ggtitle("Estimated times of strong speciation rate shifts - RPANDA models")+
  theme(legend.position = "top", axis.text.x=element_text(angle=45,hjust=1),plot.title = element_text(hjust = 0.5))


d_shift_rpanda <- sim.d.shift.matrix %>% pivot_longer(cols = -V101,values_to = "Time") %>% dplyr::rename(Lineage=V101) %>% 
  ggplot(aes(x=Lineage,y=Time))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  ylim(c(0.5,160))+
  theme_bw()+
  ggtitle("Estimated times of strong extinction rate shifts - RPANDA models")+
  theme(legend.position = "top", axis.text.x=element_text(angle=45,hjust=1),plot.title = element_text(hjust = 0.5))

b_shift_comet <- ggplot(comet.b.shift.times,aes(x=Lineage,y=Time))+
  #geom_boxplot(outliers = F)+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_bw()+
  ggtitle("Estimated times of strong speciation rate shifts - CoMET")+
  theme(legend.position = "top", axis.text.x=element_text(angle=45,hjust=1),plot.title = element_text(hjust = 0.5))+
  geom_point(data=data, aes(y=Up_Peak_or_Down_Dip_Time),fill="blue",shape=24,color="white",size=3.5,alpha=0.7)+
  geom_point(data=data, aes(y=Peak_or_Dip_Time),fill="blue",shape=21,color="white",size=3.5,alpha=0.7)+
  geom_point(data=data, aes(y=Down_Peak_or_Up_Dip_Time),fill="blue",shape=25,color="white",size=3.5,alpha=0.7)



grid.arrange(b_shift_rpanda,d_shift_rpanda,b_shift_comet,nrow=2)
dev.off()

# set the working directory back to the root of the R project
setwd("../../../")