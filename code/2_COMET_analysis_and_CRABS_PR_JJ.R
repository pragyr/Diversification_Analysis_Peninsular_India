library(ape)
library(TESS)
library(tibble)
library(pracma)
library(dplyr)
library(ggplot2)
library(CRABS)


# MAIN CODE STARTS HERE
# Provide all your tree files in a folder with the name - 'all_trees' and then set the working directory to the path of that 'all_trees' folder


setwd("data/all_trees")

system2("mkdir",args=c("../COMET_results"))
system2("mkdir",args=c("../TESS_RTT_data"))
system2("mkdir",args=c("../CRABS_results/"))

filenames <- list.files()
lineage.names <- substr(filenames,1,(nchar(filenames)-4))
mean.rates <- matrix(ncol=3)
colnames(mean.rates) <- c("Lineage names",
                          "Speciation rates",
                          "Extinction rates")

for(i in 1:length(filenames)){
  #i <- 2
  myTree <- read.tree(filenames[i])
  if(!is.ultrametric(myTree)) {myTree <- force.ultrametric(myTree,"extend")}
  cat("\nCurrently using tree of",lineage.names[i],"\n")
  times <- as.numeric(branching.times(myTree))
  cr.age <- max(times)
  ##### rjMCMC and CoMET analysis with empirical hyperpriors
  numExpectedMassExtinctions <- 1
  numExpectedRateChanges <- 1
  
  samplingFraction <- 1
  
  # Make directory with same name as tree
  system2("mkdir",args=c(paste0("../COMET_results/",lineage.names[i])))
  setwd(paste0("../COMET_results/",lineage.names[i]))
  
  # Perform TESS-COMET analysis
  convergence <- FALSE
  niter <- 1000
  while(!convergence){
    tess.analysis(myTree,
                  empiricalHyperPriors = TRUE,
                  samplingProbability = samplingFraction,
                  estimateNumberMassExtinctions = FALSE,
                  BURNIN = (niter/10),
                  MAX_ITERATIONS = niter,
                  dir = "COMET_without_mass_ex")
    
    print(c(paste0(lineage.names[i],"_comet_no_mass_extinctions done")))
    
    output <- tess.process.output("COMET_without_mass_ex",
                                  numExpectedRateChanges = numExpectedRateChanges,
                                  numExpectedMassExtinctions = numExpectedMassExtinctions)
    sp.shift.df <- as.data.frame(output$`speciation shift times`)
    ex.shift.df <- as.data.frame(output$`extinction shift times`)
    
    bayes.factors <- as.data.frame(rbind(output$`speciation Bayes factors`,output$`extinction Bayes factors`))
    
    sp.rates.df <- as.data.frame(output$`speciation rates`)
    ex.rates.df <- as.data.frame(output$`extinction rates`)
    colnames(sp.rates.df) <- colnames(ex.rates.df) <- colnames(sp.shift.df) <- colnames(ex.shift.df) <- colnames(bayes.factors) <- output$intervals[-length(output$intervals)]
    
    
    ess.vals <- c(
      effectiveSize(output$numSpeciationCategories),
      effectiveSize(output$numExtinctionCategories))
    
    if(length(which(ess.vals <= 200)) == 0){
      convergence <- TRUE
      cat("\nConverged. Writing the results.\n\n")
      mean.sp.rate.across.time <- mean(unlist(lapply(sp.rates.df, median)))
      mean.ex.rate.across.time <- mean(unlist(lapply(ex.rates.df, median)))
      mean.rates <- rbind(mean.rates,c(lineage.names[i],mean.sp.rate.across.time,mean.ex.rate.across.time))
      
      write.csv(sp.rates.df,file=paste0("../../TESS_RTT_data/",lineage.names[i],"_sp_rates.csv"))
      write.csv(ex.rates.df,file=paste0("../../TESS_RTT_data/",lineage.names[i],"_ex_rates.csv"))
      write.csv(sp.shift.df,file=paste0("../../TESS_RTT_data/",lineage.names[i],"_sp_shift_times.csv"))
      write.csv(ex.shift.df,file=paste0("../../TESS_RTT_data/",lineage.names[i],"_ex_shift_times.csv"))
      write.csv(bayes.factors,file=paste0("../../TESS_RTT_data/",lineage.names[i],"_bayes_factors.csv"))
    }else{
      cat("\nDid not converge. Re-running with 5 times the previous number of iterations.\n")
      niter = niter * 5
    }
  }
  
  ### PLOTS TO CHECK FOR CONVERGENCE
  pdf(paste0(lineage.names[i],"_comet_no_mass_extinctions_convergence.pdf"), 6, 6, bg="transparent")
  layout.mat <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
  layout(layout.mat)
  tess.plot.singlechain.diagnostics(output,
                                    parameters = c("speciation rates",
                                                   "extinction rates"),
                                    las=2)
  
  dev.off()
  
  
  
  ### PLOTS FOR RATES
  pdf(paste0(lineage.names[i],"_speciation_and_extinction_rates.pdf"), 6, 6, bg="transparent")
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
  
  mtext(lineage.names[i], outer=TRUE, cex = 2.0, col="Black")
  
  dev.off()
  
  ############# CRABS TEST ################
  est_speciation <- apply(output[["speciation rates"]], 2, median)
  est_extinction <- apply(output[["extinction rates"]], 2, median)
  
  times <- seq(0, cr.age, length.out = 100)
  
  ebd_tess <- tibble(
    "time" = times,
    "lambda" = rev(est_speciation),
    "mu" = rev(est_extinction),
  )
  
  
  ############ The following snippet of code has been taken and modified from Vivek Cyriac's paper -https://doi.org/10.1093/evolut/qpae006 ##################
  ## CRABS models
  ## settings
  ## Calculate maximum/min slope
  slopem <- function(f, times, method = "max"){
    fd <- fderiv(Vectorize(f), times, method = "central")
    res <- eval(call(method, fd))
    return(res)
  }
  
  logistic_shift <- function(min, max, steepness, x0){
    res <- function(x) {
      ((max-min) / (1 + steepness^(-1*(x - x0)))) + min
    }
    return(res)
  }
  
  lseq <- function (from = 1, to = 1e+05, length.out = 6) {
    exp(seq(log(from), log(to), length.out = length.out))
  }
  
  expseq <- function (from = 1, to = 1e+05, length.out = 6) {
    log(seq(exp(from), exp(to), length.out = length.out))
  }
  
  concat.factor <- function (...) {
    as.factor(do.call(c, lapply(list(...), as.character)))
  }
  
  height <- max(node.depth.edgelength(myTree))
  times <- seq(0, height, length.out = 500)
  n <- 5
  
  foo_rate <- list()
  
  ## Constants
  bs <- seq(0.1, 0.6, length.out = n)
  constant <- lapply(bs, function(b) function(t) b)
  names(constant) <- paste0("constant", seq_along(constant))
  foo_rate[["constant"]] <- constant
  
  
  ## Linear-up
  bs <- seq(0.001, 0.007, length.out = n)
  linear_up <- sapply(bs, function(b) function(t) 0.6 - b * t)
  names(linear_up) <- paste0("linear_up", seq_along(linear_up))
  foo_rate[["linear_up"]] <- linear_up
  
  ## Linear-down
  bs <- seq(0.001, 0.009, length.out = n)
  linear_down <- sapply(bs, function(b) function(t) 0.1 + b * t)
  names(linear_down) <- paste0("linear_down", seq_along(linear_down))
  foo_rate[["linear_down"]] <- linear_down
  
  ## Exponentially increasing
  bs <- seq(0.05, 0.2, length.out = n)
  exp_up <- sapply(bs, function(b) function(t) 0.1 + 0.3*exp(-b*t))
  names(exp_up) <- paste0("exp_up", seq_along(exp_up))
  foo_rate[["exp_up"]] <- exp_up
  
  ## Exponential decreasing
  exp_down <- sapply(bs, function(b) function(t) 0.1 + 0.3*exp(-b*(height-t)))
  names(exp_down) <- paste0("exp_down", seq_along(exp_down))
  foo_rate[["exp_down"]] <- exp_down
  
  ## Sigmoidal Up-shift
  steepnesses <- lseq(0.1, 0.9, length.out = n)
  up <- lapply(steepnesses, function(s) logistic_shift(min=0.1, max=0.6, steepness = s, x0 = 2.5))
  names(up) <- paste0("up", seq_along(up))
  foo_rate[["up"]] <- up
  
  sapply(foo_rate[["up"]], function(f) slopem(f, times, method = "min"))
  
  ## Signmoidal Down-shift
  #steepnesses2 <- lseq(1.12, 5, length.out = n)
  mse_slope <- function(f, times, target){
    s <- slopem(f, times, method = "max")
    mse <- (target - s)^2
    return(mse)
  }
  
  steepnesses2l <- list()
  for (j in 1:n){
    slope_up <- (-1)*slopem(foo_rate[["up"]][[j]], times, method = "min")
    steepnesses2l[[j]] <- optimize(function(s) mse_slope(logistic_shift(min=0.1, max=0.6, steepness = s, x0 = 2.5), times, slope_up), c(1.01, 15.0))$minimum
  }
  steepnesses2 <- unlist(steepnesses2l)
  
  down <- lapply(steepnesses2, function(s) logistic_shift(min=0.1, max=0.6, steepness = s, x0 = height/2.0))
  names(down) <- paste0("down", seq_along(down))
  foo_rate[["down"]] <- down
  
  ## now the maximum slope of the S- rates are equal to the minimum slopes of the S+ (times *(-1)), in reversed order
  rev(sapply(foo_rate[["down"]], function(f) slopem(f, times, method = "max")))
  
  
  l <- list()
  for (j in seq_along(foo_rate)){
    item <- names(foo_rate)[j]
    print(item)
    foos <- foo_rate[[j]]
    
    df <- as_tibble(sapply(foos, function(foo) sapply(times, foo)))
    df$time <- times
    
    df <- gather(df, key = "subitem", value = "rate", -time)
    df$item <- item
    df$subitemidx <- readr::parse_number(df$subitem)
    df$is_three <- df$subitemidx == 3
    
    l[[j]] <- df
  }
  pdata <- bind_rows(l)
  
  p <- ggplot(pdata, aes(x = time, y = rate, linetype = factor(subitemidx), color = factor(is_three))) +
    facet_wrap(.~ item, nrow = 2) +
    geom_line() +
    theme_classic() +
    scale_x_reverse() +
    scale_color_manual(values = c("gray", "blue")) +
    scale_linetype_manual(values = rep("solid", 6))
  
  
  ## Reference model
  ## What is the average extinction rate of the proposed rate functions?
  avg_lambda <- sapply(unlist(foo_rate), function(f) quadgk(f, 0, height)/height)
  
  median(avg_lambda)
  mean(avg_lambda)
  
  ref_foos <- list(
    "constant" = foo_rate$constant$constant3,
    "linear_up" = foo_rate$linear_up$linear_up3,
    "linear_down" = foo_rate$linear_down$linear_down3,
    "up" = foo_rate$up$up3,
    "down" = foo_rate$down$down3,
    "exp_up" = foo_rate$exp_up$exp_up3,
    "exp_down" = foo_rate$exp_down$exp_down3
  )
  
  lambda <- approxfun(ebd_tess$time, ebd_tess$lambda)
  mu <-approxfun(ebd_tess$time, ebd_tess$mu)
  max_t <- max(ebd_tess[["time"]])
  times_fine <- seq(0, max_t, length.out = 500)
  my_model <- create.model(lambda,  mu, times_fine)
  
  
  references <- lapply(ref_foos, function(foos) create.model(func_spec0 = lambda, func_ext0 = foos, times = times))
  
  ## CONGRUENCE CLASSES
  cgs_fromlambda <- lapply(references, function(ref) congruent.models(my_model, mus = unlist(foo_rate), keep_ref = TRUE))
  cgs <- cgs_fromlambda
  group_names <- factor(c(names(foo_rate), "reference"), levels = c(names(foo_rate), "reference"))
  ps_lambdas <- lapply(cgs,
                       function(cg) summarize.trends(cg, threshold = 0.02, 
                                                     rate_name = "lambda",
                                                     group_names = group_names)); 
  ps_mus <- lapply(cgs,
                   function(cg) summarize.trends(cg, threshold = 0.02, 
                                                 rate_name = "mu",
                                                 group_names = group_names))
  ps_deltas <- lapply(cgs,
                      function(cg) summarize.trends(cg, threshold = 0.02, 
                                                    rate_name = "delta",
                                                    group_names = group_names))
  
  names(ps_lambdas) <- names(ps_mus) <- names(ps_deltas) <- names(foo_rate)
  
  
  pdf(paste0("../../CRABS_results/",lineage.names[i],"_CRABS_output_summary.pdf"), 6, 6, bg="transparent")
  plot(p)
  plot(my_model)
  print(ps_lambdas)
  print(ps_mus)
  print(ps_deltas)
  dev.off()
  
  setwd("../../all_trees")
}
# get mean rates through time bins
mean.rates <- mean.rates[-1,]
mean.rates <- cbind(mean.rates,
                    "Net-diversification rates"=(as.numeric(mean.rates[,"Speciation rates"]) - 
                                                   as.numeric(mean.rates[,"Extinction rates"])))
write.csv(mean.rates,"../Mean_rates_CoMET.csv",row.names = F)

# extract bayes factors through time bins
bayes_factor_files <- list.files("../TESS_RTT_data/",pattern="bayes")
# find out times of strong rate-shifts
sp.rate.shift.times <- ex.rate.shift.times <- rep(NA, length(bayes_factor_files))
names(sp.rate.shift.times) <- names(ex.rate.shift.times) <- lineage.names

# loop through the time bins and find out the ones where bayes factor is >4.6 (i.e. support for episodic rate shift is high enough)
# check methods for justification of the threshold
for(i in 1:length(bayes_factor_files)){
  #i=23
  bf <- read.csv(paste0("../TESS_RTT_data/",bayes_factor_files[i]))[,-1]
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
write.csv(cbind("Lineage names"=lineage.names,
                "Speciation rate shift times"=as.vector(sp.rate.shift.times),
                "Extinction rate shift times"=as.vector(ex.rate.shift.times)),"../times_of_strong_rate_shifts_CoMET.csv")

######################################### TESS (CoMET) data editing - add net-div data ###############################################
tess.files <- list.files("../TESS_RTT_data/")
tess.fn <- unique(substr(tess.files,1,(nchar(tess.files)-13)))

for(i in 1:length(tess.fn)){
  #i=1
  times.dat <- colnames(read.csv(paste0("../TESS_RTT_data/",tess.fn[i],"_sp_rates.csv")))
  sp.rates <- colMeans(read.csv(paste0("../TESS_RTT_data/",tess.fn[i],"_sp_rates.csv"))[,-1])
  ex.rates <- colMeans(read.csv(paste0("../TESS_RTT_data/",tess.fn[i],"_ex_rates.csv"))[,-1])
  nd.rates <- as.matrix(t(as.vector(sp.rates) - as.vector(ex.rates)))
  
  colnames(nd.rates) <- times.dat[-1]
  write.csv(nd.rates,paste0("../TESS_RTT_data/",tess.fn[i],"_nd_rates.csv"))
}

# set the working directory back to the root of the R project
setwd("../../")