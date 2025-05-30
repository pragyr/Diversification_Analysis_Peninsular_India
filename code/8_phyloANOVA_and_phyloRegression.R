library(ape)
library(MASS)
library(mvtnorm)
library(caper)
library(RRPP)
library(multcomp)
library(nlme)
library(maps)
library(phytools)
library(FSA)

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(olsrr)

library(phyr)


library(geiger)



library(ggstatsplot)
library(tidytree)

library(ggdist)

setwd("data/")

options(scipen = 100,digits=7)


data1 <- read.csv("Final_summary_updated.csv")

tree1 = read.tree("SUPER_TREE.nwk")

# The non-parametric Kruskal-Wallis Test and the Dunn Test
kruskal_wallis_and_dunntest <- function(formula,dat,p.adj.meth="bonferroni",comp.name){
  #dat <- data;formula <- Net.diversification.rate..CoMET. ~  Taxonomic.group;comp.name = "Diversification rate vs. Taxonomic group"
  m.kruskal <- kruskal.test(formula,data = dat)
  
  res.dunntest <- dunnTest(formula,data = dat, method = p.adj.meth)
  p.vals <- c(m.kruskal$p.value,unlist(unname(res.dunntest$res[c("P.adj")])))
  names(p.vals) <- c(paste("Overall Kruskal-Wallis p-value for",comp.name),paste(unlist(unname(res.dunntest$res["Comparison"]))))
  return(p.vals)
}

# Phylogenetic Generalised Least-Squares Regression
pgls_phyloanova <- function(formula,is.predictor.categorical=FALSE,phy,df.pgls){
  #formula=Net.diversification.rate..RPANDA.~Taxonomic.group;is.predictor.categorical=T;phy=tree1;df.pgls=data1
  Response <- as.character(formula)[2]
  Predictor <- as.character(formula)[3]
  
  colnames(df.pgls)[which(colnames(df.pgls) == Response)] <- "Response"
  colnames(df.pgls)[which(colnames(df.pgls) == Predictor)] <- "Predictor"
  
  
  df.pgls <- df.pgls[match(phy$tip.label,if("tip_names" %in% colnames(df.pgls)){df.pgls$tip_names}else{df.pgls$Lineage}),] # order data in the same format as tip labels of tree
  
  comp.data <- if("tip_names" %in% colnames(df.pgls)){
    comparative.data(phy=phy,data = df.pgls,names.col = tip_names,vcv = TRUE, 
                     na.omit = FALSE, warn.dropped = TRUE)
  }else{
    comparative.data(phy=phy,data = df.pgls,names.col = Lineage,vcv = TRUE, 
                     na.omit = FALSE, warn.dropped = TRUE)
  }
  
  
  pgls.model <- pgls(Response~Predictor,data = comp.data,lambda = 'ML')
  
  summ.pgls <- summary(pgls.model) # summary of the pgls model
  lambda.r <- summ.pgls$param["lambda"] # value for the lambda parameter
  lambda.r.p <- summ.pgls$param.CI$lambda$bounds.p[1] # p value for the value of lambda
  coeff.r <- as.data.frame(summ.pgls$coefficients) # coefficients of the pgls model
  p.val.r <- coeff.r$`Pr(>|t|)`[-1] # p-value of the model (i.e. for the slope)
  names(p.val.r) <- rownames(coeff.r)[-1]
  intercept.r <- coeff.r[1,"Estimate"]
  slope.r <- coeff.r[2,"Estimate"]
  
  
  r2.r <- summ.pgls$r.squared # multiple r squared value
  
  stats <- c("R squared"=r2.r,
             "p-value"=round(unname(p.val.r[1]),5),
             "lambda"=lambda.r,
             "lambda p-value"=lambda.r.p)
  

  if(is.predictor.categorical){
    lineage.or.tip.names <- if("tip_names" %in% colnames(df.pgls)){df.pgls$tip_names}else{df.pgls$Lineage}
    rownames(df.pgls) <- lineage.or.tip.names
    
    df.pgls$Predictor.mod <- as.factor(df.pgls$Predictor)
    Predictor.mod <- df.pgls$Predictor.mod
    names(Predictor.mod) <- lineage.or.tip.names
    Response <- df.pgls$Response
    names(Response) <- lineage.or.tip.names
    # make correlation structure from phylogeny
    corBM<-corBrownian(phy=phy,form=~lineage.or.tip.names)
    cov.mat <- vcv(corBM)
    
    cat("\nPerforming Phylogenetic ANOVA\n")
    phylo.anova <- lm.rrpp(Response~Predictor.mod,
                           #data = rrpp.data.frame(df.pgls),
                           Cov = cov.mat,
                           print.progress = T,
                           verbose = T)
    
    
    sum.phylo.anova <- summary(phylo.anova)
    phylo.anova.table <- sum.phylo.anova$table[c("Df","SS","Residual SS","Rsq","F","Pr(>F)")]
    names(phylo.anova.table) <- c("df","SSR","SSE","MSR","F-value","p-value")
    
    cat("\nPerforming Post-Hoc cmoparisons\n")
    phyloANOVA.post.hoc <- summary(pairwise(phylo.anova,groups = Predictor.mod))$summary.table
    names <- rownames(phyloANOVA.post.hoc)
    phyloANOVA.post.hoc <- phyloANOVA.post.hoc$`Pr > d`
    names(phyloANOVA.post.hoc) <- gsub(":"," - ",names)
    
    
    return(list("stats"=stats,"lambda"=lambda,
                "phyloANOVA.results.phylANOVA"=phylo.anova.table,
                "phyloANOVA.post.hoc"=phyloANOVA.post.hoc))
  }else{
    return(list("stats"=stats))
  }
  
}

# Phylogenetic signal in residuals
phylosig_in_residuals <- function(formula,dat,phy){
  lm.model <- lm(formula = formula,data=dat)
  residuals <- residuals(lm.model)
  names(residuals) <- dat$Lineage
  residuals <- residuals[match(phy$tip.label,names(residuals))]
  lambda.res <- phylosig(tree = phy,x = residuals,method = "lambda",test = T)
  return(list("lambda"=lambda.res$lambda,"lambda.p.val"=lambda.res$P))
}

summary_info <- function(response,data,tree){
  stats.rate.vs.tg <- list(kruskal_wallis_and_dunntest(as.formula(paste(response,"~","Taxonomic.group")),data,comp.name = "Diversification rate vs. Taxonomic group"),
                           pgls_phyloanova(as.formula(paste(response,"~","Taxonomic.group")),T,tree,data),
                           phylosig_in_residuals(as.formula(paste(response,"~","Taxonomic.group")),data,tree))
  
  stats.rate.vs.bo <- list(kruskal_wallis_and_dunntest(as.formula(paste(response,"~","Biogeographic.origin")),data,comp.name = "Diversification rate vs. Biogeographic origin"),
                           pgls_phyloanova(as.formula(paste(response,"~","Biogeographic.origin")),T,tree,data),
                           phylosig_in_residuals(as.formula(paste(response,"~","Biogeographic.origin")),data,tree))
  
  stats.rate.vs.ht <- list(kruskal_wallis_and_dunntest(as.formula(paste(response,"~","Habitat")),data,comp.name = "Diversification rate vs. Habitat type"),
                           pgls_phyloanova(as.formula(paste(response,"~","Habitat")),T,tree,data),
                           phylosig_in_residuals(as.formula(paste(response,"~","Habitat")),data,tree))
  
  if(!("tip_names" %in% colnames(data))){
    stats.rate.vs.ca <- list(summary(lm(as.formula(paste(response,"~","Crown.age..Mya.")),data)),
                             pgls_phyloanova(as.formula(paste(response,"~","Crown.age..Mya.")),F,tree,data),
                             phylosig_in_residuals(as.formula(paste(response,"~","Crown.age..Mya.")),data,tree))
  }
  
  phyloanova_pgls.stats <- rbind("Taxonomic groups (PhyloANOVA)"=stats.rate.vs.tg[[2]]$phyloANOVA.results.phylANOVA,
                            "Biogeographic origins (PhyloANOVA)"=stats.rate.vs.bo[[2]]$phyloANOVA.results.phylANOVA,
                            "Habitat types (PhyloANOVA)"=stats.rate.vs.ht[[2]]$phyloANOVA.results.phylANOVA)
  if(!("tip_names" %in% colnames(data))){
    phyloanova_pgls.stats <- list(phyloanova_pgls.stats,
                                  "Clade age (PGLS)"=stats.rate.vs.ca[[2]]$stats)
  }
  
  kruskal_lm.stats <- rbind("Taxonomic groups (KW)"=stats.rate.vs.tg[[1]][1],
                       "Biogeographic origins (KW)"=stats.rate.vs.bo[[1]][1],
                       "Habitat types (KW)"=stats.rate.vs.ht[[1]][1])
  if(!("tip_names" %in% colnames(data))){
    kruskal_lm.stats <- list(kruskal_lm.stats,
                                  "Clade age (LM)"=c("R squared"=stats.rate.vs.ca[[1]]$adj.r.squared,
                                    "p-value"=unname(stats.rate.vs.ca[[1]]$coefficients[,"Pr(>|t|)"][2])))
  }
  
  phyloanova_pgls.post.hoc.p.values <- c(stats.rate.vs.tg[[2]]$phyloANOVA.post.hoc,
                           stats.rate.vs.ht[[2]]$phyloANOVA.post.hoc)
  
  
  kruskal_lm.post.hoc.p.values <- c(stats.rate.vs.tg[[1]][c(2,3,4)],
                      stats.rate.vs.ht[[1]][c(2,3,4)])
  
  phylogentic.signal <- rbind("Residuals\n(Rate vs. Taxonomic group)"=unlist(stats.rate.vs.tg[[3]][c("lambda","lambda.p.val")]),
                              "Residuals\n(Rate vs. Biogeographic origin)"=unlist(stats.rate.vs.bo[[3]][c("lambda","lambda.p.val")]),
                              "Residuals\n(Rate vs. Habitat type)"=unlist(stats.rate.vs.ht[[3]][c("lambda","lambda.p.val")]))
  return(list("PhyloANOVA stats"=phyloanova_pgls.stats,"PhyloANOVA post hoc p values"=phyloanova_pgls.post.hoc.p.values,
              "KW stats"=kruskal_lm.stats,"KW post hoc (Dunn's Test) p values"=kruskal_lm.post.hoc.p.values,"Phylogenetic signal"=phylogentic.signal))
  
}
summary.RPANDA <- summary_info(response="Net.diversification.rate..RPANDA.",data1,tree1)
summary.CoMET <-  summary_info(response="Net.diversification.rate..CoMET.",data1,tree1)
summary.ClaDS <-  summary_info(response="Net.diversification.rate..ClaDS.",data1,tree1)


###################### TIP RATES ####################################
data2 <- read.csv("clads_tip_rates_with_attributes.csv")
tree2 <- read.tree("SUPER_TREE_with_all_tips_ultrametric.nwk")
summary_info(response = "Net.diversification_Rate",data = data2,tree=tree2)

# set the working directory back to the root of the R project
setwd("../")