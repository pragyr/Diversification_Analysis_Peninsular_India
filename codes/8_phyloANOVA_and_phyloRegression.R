library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(olsrr)
library(FSA)
library(phyr)
library(ape)
library(phytools)
library(geiger)
library(nlme)
library(multcomp)
library(caper)
library(gridExtra)
library(ggstatsplot)
library(tidytree)
library(nortest)
library(RRPP)
library(ggdist)
library(adephylo)
library(phylobase)
library(phylosignal)

setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../data/"))

options(scipen = 100,digits=7)


data <- read.csv("Final_summary_updated2.csv")
phylo_DIV <- function(name){
  phy <- read.tree(paste0("./all_trees/",name))
  pd <- sum(phy$edge.length)
  names(pd) <- substr(name,1,nchar(name)-4)
  return(pd)
}
pds <- unlist(lapply(list.files("all_trees"), phylo_DIV))
pds <- pds[-which(names(pds)=="Montecincla")]
length(pds)
length(data$No..of.species)
m <- lm(data$No..of.species ~ pds)

data$Phylogenetic_diversity <- pds
write.csv(data,"Final_summary_updated2.csv")
summary(m)
plot(m)
tree = read.tree("SUPER_TREE.nwk")
tree = drop.tip(tree,"Montecincla")

kruskal_wallis_and_dunntest <- function(formula,dat,p.adj.meth="bonferroni",comp.name){
  #dat <- data;formula <- Net.diversification.rate..CoMET. ~  Taxonomic.group;comp.name = "Diversification rate vs. Taxonomic group"
  m.kruskal <- kruskal.test(formula,data = dat)
  res.dunntest <- dunnTest(formula,data = dat, method = p.adj.meth)
  p.vals <- c(m.kruskal$p.value,unlist(unname(res.dunntest$res[c("P.adj")])))
  names(p.vals) <- c(paste("Overall Kruskal-Wallis p-value for",comp.name),paste("Dunn test p-value for",unlist(unname(res.dunntest$res["Comparison"]))))
  return(p.vals)
}

anova_with_post_hoc <- function(formula,dat,p.adj.meth="bonferroni",comp.name){
  #dat <- data;formula <- Net.diversification.rate..CoMET. ~  Taxonomic.group;comp.name = "Diversification rate vs. Taxonomic group"
  m.aov <- aov(formula = formula,data = dat)
  sum.aov <- as.data.frame(summary(m.aov)[[1]])
  sum.aov <- unlist(sum.aov[c("Df","Sum Sq","Mean Sq","F value","Pr(>F)")])
  sum.aov <- sum.aov[which(!is.na(sum.aov))]
  sum.aov <- sum.aov[-c(2,6)]
  
  names(sum.aov) <- c("df","SSR","SSE","MSR","F-value","p-value")
  
  res.tukeyhsd <- TukeyHSD(m.aov)
  
  p.vals <- res.tukeyhsd[[1]][,"p adj"]
  names(p.vals) <- paste("TukeyHSD test p-value for",names(p.vals))
  return(list("ANOVA_stats"=sum.aov,"TukeyHSD_pairiwise_comparisons"=p.vals))
}

convert.diag.matrix.to.vector <- function(mat){
  vec <- c()
  name.vec <- c()
  for(i in 1:nrow(mat)){
    for(j in 1:i){
      if(j==i) next
      vec <- c(vec,mat[i,j])
      name.vec <- c(name.vec,paste0(rownames(mat)[i],"-",colnames(mat)[j]))
    }
  }
  names(vec) <- name.vec
  return(vec)
}

mk_pgls_plot <- function(formula,is.predictor.categorical=FALSE,phy,df.pgls,xlab,ylab,make.plot=TRUE,xlim.min,ylim.min){
  # formula = Net.diversification.rate..CoMET. ~  Taxonomic.group
  # is.predictor.categorical=T;phy=tree;df.pgls = data;make.plot = F
  # xlim.min=0;ylim.min=0
  # formula <- Net.diversification.rate..CoMET. ~ Habitat
  # is.predictor.categorical <- T
  # phy <- tree
  # df.pgls <- data
  # make.plot = FALSE
  formula=Net.diversification.rate..ClaDS.~Crown.age..Mya.
  is.predictor.categorical=FALSE
  phy=tree
  df.pgls=data
  
  
  Response <- as.character(formula)[2]
  Predictor <- as.character(formula)[3]
  
  colnames(df.pgls)[which(colnames(df.pgls) == Response)] <- "Response"
  colnames(df.pgls)[which(colnames(df.pgls) == Predictor)] <- "Predictor"
  
  
  df.pgls <- df.pgls[match(phy$tip.label,df.pgls$Lineage),] # order data in the same format as tip labels of tree
  
  comp.data <- comparative.data(phy=phy,data = df.pgls,names.col = Lineage,vcv = TRUE, 
                                na.omit = FALSE, warn.dropped = TRUE)
  
  lm.model <- lm(Response~Predictor,data=df.pgls)
  residuals.lm <- residuals(lm.model)
  names(residuals.lm) <- df.pgls$Lineage
  residuals.lm <- residuals.lm[match(phy$tip.label,names(residuals.lm))]
  lambda.res.lm <- phylosig(tree = phy,x = residuals.lm,method = "lambda",test = T)
  
  summary(lm.model)
  lambda.res.lm
  
  
  
  
  pgls.model <- pgls(Response~Predictor,data = comp.data,lambda = 'ML')
  residuals.pgls <- residuals(pgls.model)
  names.res <- rownames(residuals.pgls)
  residuals.pgls <- c(residuals.pgls)
  names(residuals.pgls) <- names.res
  residuals.pgls <- residuals.pgls[match(phy$tip.label,names(residuals.pgls))]
  lambda.res <- phylosig(tree = phy,x = residuals.pgls,method = "lambda",test = T)
  
  summary(pgls.model)
  lambda.res
  
  
  # pgls.model1 <- pgls(Net.diversification.rate..ClaDS. ~ Crown.age..Mya.+Taxonomic.group+Biogeographic.origin+Habitat,
  #                     data=comp.data,lambda = 'ML')
  # pgls.model2 <- pgls(Net.diversification.rate..ClaDS. ~ Crown.age..Mya.+Taxonomic.group+Biogeographic.origin,
  #                     data=comp.data,lambda = 'ML')
  # pgls.model3 <- pgls(Net.diversification.rate..ClaDS. ~ Crown.age..Mya.+Taxonomic.group+Habitat,
  #                     data=comp.data,lambda = 'ML')
  # pgls.model4 <- pgls(Net.diversification.rate..ClaDS. ~ Crown.age..Mya.+Taxonomic.group,
  #                     data=comp.data,lambda = 'ML')
  # pgls.model5 <- pgls(Net.diversification.rate..ClaDS. ~ Crown.age..Mya.+Biogeographic.origin+Habitat,
  #                     data=comp.data,lambda = 'ML')
  # pgls.model6 <- pgls(Net.diversification.rate..ClaDS. ~ Crown.age..Mya.+Biogeographic.origin,
  #                     data=comp.data,lambda = 'ML')
  # pgls.model7 <- pgls(Net.diversification.rate..ClaDS. ~ Crown.age..Mya.+Habitat,
  #                     data=comp.data,lambda = 'ML')
  # pgls.model8 <- pgls(Net.diversification.rate..ClaDS. ~ Crown.age..Mya.,
  #                     data=comp.data,lambda = 'ML')
  # 
  # summ1 <- summary(pgls.model1)
  # summ2 <- summary(pgls.model2)
  # summ3 <- summary(pgls.model3)
  # summ4 <- summary(pgls.model4)
  # summ5 <- summary(pgls.model5)
  # summ6 <- summary(pgls.model6)
  # summ7 <- summary(pgls.model7)
  # summ8 <- summary(pgls.model8)
  # 
  # c(summ1$r.squared,summ2$r.squared,summ3$r.squared,
  #   summ4$r.squared,summ5$r.squared,summ6$r.squared,
  #   summ7$r.squared,summ8$r.squared)
  # c(summ1$coefficients)
  
  summ.pgls <- summary(pgls.model)
  lambda.r <- summ.pgls$param["lambda"]
  coeff.r <- as.data.frame(summ.pgls$coefficients)
  p.val.r <- coeff.r$`Pr(>|t|)`[-1]
  names(p.val.r) <- rownames(coeff.r)[-1]
  intercept.r <- coeff.r[1,"Estimate"]
  slope.r <- coeff.r[2,"Estimate"]
  
  r2.r <- summ.pgls$r.squared
  stats <- as.character(as.expression(substitute(italic("R")^2~"="~r2*","~italic("p")~"="~p,
                                                 list(r2 = round(r2.r,4),
                                                      p = round(unname(p.val.r[1]),5)
                                                 ))))
  if(is.na(lambda.r)){
    lambda <- as.character(as.expression(substitute(italic("\u03BB")~"="~l,
                                                    list(l = "NA"))))
  }else{
    lambda <- as.character(as.expression(substitute(italic("\u03BB")~"="~l,
                                                    list(l = round(unname(lambda.r)[1],3)))))
  }
  
  
  # eq <- as.character(as.expression(substitute(italic(y) == a + b %.% italic("* x"), 
  #                                             list(a = format(intercept.r, digits = 2),
  #                                                  b = format(slope.r, digits = 2)))))
  
  
  
  if(!make.plot){
    if(is.predictor.categorical){
      lineage.names <- df.pgls$Lineage
      rownames(df.pgls) <- lineage.names
      
      df.pgls$Predictor.mod <- as.factor(df.pgls$Predictor)
      Predictor.mod <- df.pgls$Predictor.mod
      names(Predictor.mod) <- lineage.names
      Response <- df.pgls$Response
      names(Response) <- lineage.names
      # make variance covariance matrix from phylogeny
      cov.mat <- vcv.phylo(phy = tree)
      cat("\nPerforming GLS\n")
      phylo.anova <- lm.rrpp(Response~Predictor.mod,
              #data = rrpp.data.frame(df.pgls),
              Cov = cov.mat,
              print.progress = T,
              verbose = T)
      
      sum.phylo.anova <- summary(phylo.anova)
      phylo.anova.table <- sum.phylo.anova$table[c("Df","SS","Residual SS","Rsq","F","Pr(>F)")]
      names(phylo.anova.table) <- c("df","SSR","SSE","MSR","F-value","p-value")
      
      
      if(length(unique(Predictor.mod))>2){
        cat("\nPerforming Post-Hoc\n")
        summary.pairwise.comp <- summary(pairwise(fit = phylo.anova,groups = Predictor.mod,
                                                  print.progress = T,verbose = T))
        
        phylANOVA.post.hoc <- convert.diag.matrix.to.vector(summary.pairwise.comp$pairwise.tables$P)
      }else{
        phylANOVA.post.hoc <- NA
      }
      
      
      return(list("stats"=stats,"lambda"=lambda,
                  "phyloANOVA.results.phylANOVA"=phylo.anova.table,
                  "phyloANOVA.post.hoc"=phylANOVA.post.hoc))
    }else{
      return(list("stats"=stats,"lambda"=lambda))
    }
    
  }else{
    if(p.val.r > 0.05){
      lty<-3
    }else{lty<-1}
    colnames(df.pgls)[which(colnames(df.pgls) == Response)] <- "Response"
    colnames(df.pgls)[which(colnames(df.pgls) == Predictor)] <- "Predictor"
    
    get_labels <- function(x,lab) {
      if(!(-1 %in% unlist(gregexpr("log",lab)))){
        exp(x)
      }else{
        x
      }
      
    }
    
    get_breaks_custom <- function(vals,lab){
      if(!(-1 %in% unlist(gregexpr("log",lab)))){
        log(c(1,2,5,10,20,50,100,200))
      }else{
        max.val <- max(vals)
        min.val <- min(vals)
        interval <- max.val/5
        if(min.val > 1){
          max.val <- round(max.val/5,0)*5
          interval <- round(interval,0)
        }else{
          max.val <- round(max.val/5,1)*5
          interval <- round(interval,1)
        }
        seq(0,max.val,interval)
      }
    }
    
    breaks.x <- get_breaks_custom(df.pgls$Predictor,xlab)
    breaks.y <- get_breaks_custom(df.pgls$Response,ylab)
    plot <- ggplot(df.pgls,aes(x=Predictor,y=Response)) + 
      geom_point(aes(color=Biogeographic.origin,
                     fill=Taxonomic.group,
                     shape=Habitat),size=4,stroke=0.7) +
      scale_fill_manual(values = c(alpha("#77dd77",0.8),alpha("#898989",0.8),alpha("#00c7e7",0.8)),
                         guide=guide_legend(override.aes = list(shape = 21,stroke=0,
                                                                color=c(alpha("#77dd77",0.8),alpha("#898989",0.8),alpha("#00c7e7",0.8)))))+
      scale_shape_manual(values = c(21,22,24))+
      scale_color_manual(
        values = c("#aa4231", "#222233"),
        guide = guide_legend(override.aes = list(shape = 21))
      )+
      scale_x_continuous(labels = get_labels(breaks.x,lab = xlab), breaks = breaks.x,limits = c(xlim.min,(max(df.pgls$Predictor)*1.1)))+
      scale_y_continuous(labels = get_labels(breaks.y,lab = ylab), breaks = breaks.y,limits = c(ylim.min,(max(df.pgls$Response)*1.3)))+
      xlab(ifelse(!(-1 %in% unlist(gregexpr("log",xlab))),substr(xlab,(unlist(gregexpr("log",xlab))+4),(nchar(xlab)-1)),xlab))+
      ylab(ifelse(!(-1 %in% unlist(gregexpr("log",ylab))),substr(ylab,(unlist(gregexpr("log",ylab))+4),(nchar(ylab)-1)),ylab))+
      geom_abline(intercept = intercept.r, slope = slope.r,linetype=lty) + 
      # geom_text(data = as.data.frame(eq), aes(label = eq,x = max(df.pgls$Predictor)/1.75,
      #                                         y = (max(df.pgls$Response))*1.45,family = 'serif'), parse = TRUE,size=3.5)+
      geom_text(data = as.data.frame(stats), aes(label = stats,x = max(df.pgls$Predictor)/1.1,
                                                 y = (max(df.pgls$Response))/1.5,family = 'serif'), inherit.aes = FALSE, parse = TRUE,size=3.5)+
      geom_text(data = as.data.frame(lambda), aes(label = lambda,x = max(df.pgls$Predictor)/1.2,
                                                  y = (max(df.pgls$Response))/1.8,family = 'serif'), inherit.aes = FALSE, parse = TRUE,size=3.5)+
      theme_classic()+theme(legend.title = element_blank())
    return(plot)
  }
  
}

phyloANOVA_clade_rates <- function(formula,phy,df.pgls){
  # formula = Net.diversification.rate..RPANDA. ~ Habitat
  # phy = tree
  # df.pgls <- data
  
  Response <- as.character(formula)[2]
  Predictor <- as.character(formula)[3]
  
  colnames(df.pgls)[which(colnames(df.pgls) == Response)] <- "Response"
  colnames(df.pgls)[which(colnames(df.pgls) == Predictor)] <- "Predictor"
  
  
  df.pgls <- df.pgls[match(phy$tip.label,df.pgls$Lineage),] # order data in the same format as tip labels of tree
  
  lineage.names <- df.pgls$Lineage
  rownames(df.pgls) <- lineage.names
  
  Predictor <- as.factor(df.pgls$Predictor)
  names(Predictor) <- lineage.names
  Response <- df.pgls$Response
  names(Response) <- lineage.names
  # make variance covariance matrix from phylogeny
  cov.mat <- vcv.phylo(phy = phy)
  cat("\nPerforming GLS\n")
  phylo.anova <- lm.rrpp(Response~Predictor,
                         RRPP = T,
                         Cov = cov.mat,
                         print.progress = T,
                         verbose = T)
  
  sum.phylo.anova <- summary(phylo.anova)
  phylo.anova.table <- sum.phylo.anova$table[c("Df","SS","Residual SS","Rsq","F","Pr(>F)")]
  names(phylo.anova.table) <- c("df","SSR","SSE","MSR","F-value","p-value")
  
  
  if(length(unique(Predictor.mod))>2){
    cat("\nPerforming Post-Hoc\n")
    summary.pairwise.comp <- summary(pairwise(fit = phylo.anova,groups = Predictor,
                                              print.progress = T,verbose = T))
    
    phylANOVA.post.hoc <- convert.diag.matrix.to.vector(summary.pairwise.comp$pairwise.tables$P)
  }else{
    phylANOVA.post.hoc <- NA
  }
  cat("\nDone with Post-Hoc test\n")
  residuals <- residuals(phylo.anova)
  res <- residuals[,1]
  names(res) <- rownames(residuals)
  phylosig(tree,res,method = "lambda",test=T)
  
  
  phylosig(tree,residuals(lm(Response~Predictor,df.pgls)),method = "lambda",test=T)
  
  return(list("residuals.phyloANOVA"=residuals[,1],
              "phyloANOVA.results"=phylo.anova.table,
              "phyloANOVA.post.hoc"=phylANOVA.post.hoc))
  
}

# Phylogenetic signal in residuals
phylosig_in_residuals <- function(formula,dat,phy){
  # dat <- data
  # phy <- tree
  # formula <- Net.diversification.rate..CoMET. ~ Taxonomic.group
  lm.model <- lm(formula = formula,data=dat)
  residuals <- residuals(lm.model)
  names(residuals) <- dat$Lineage
  residuals <- residuals[match(phy$tip.label,names(residuals))]
  
  lambda.res <- phylosig(tree = phy,x = residuals,method = "lambda",test = T)
  K.res <- phylosig(tree = phy,x = residuals,method = "K",test = T,nsim = 1000)
  return(list("lambda"=lambda.res$lambda,"lambda.p.val"=lambda.res$P,"K"=K.res$K,"K.p.val"=K.res$P))
}

# Phylogenetic signal in traits within groups
phylosig_in_trait_within_categories_of_a_categorical_predictor <- function(formula,dat,phy){
  # dat <- data;phy <- tree
  # formula <- Net.diversification.rate..CoMET. ~ Taxonomic.group
  Response <- as.character(formula)[2]
  Predictor <- as.character(formula)[3]
  
  categories <- unique(dat[,Predictor])
  mat.phylo.sig <- matrix(nrow=length(categories),ncol=4)
  rownames(mat.phylo.sig) <- categories
  colnames(mat.phylo.sig) <- c("lambda","lambda.p.val","K","K.p.val")
  for(g in categories){
    # g <- categories[1]
    dat.g <- dat[which(dat[,Predictor] == g),]
    phy.g <- keep.tip(phy,dat.g$Lineage)
    dat.g <- dat.g[match(phy.g$tip.label,dat.g$Lineage),]
    trait.g <- dat.g[,Response]
    names(trait.g) <- dat.g$Lineage
    mat.phylo.sig[g,c("lambda","lambda.p.val")] <- unlist(phylosig(tree = phy.g, x = trait.g,method = "lambda",test = T)[c("lambda","P")])
    mat.phylo.sig[g,c("K","K.p.val")] <- unlist(phylosig(tree = phy.g, x = trait.g,method = "K",test = T,nsim = 1000)[c("K","P")])
  }
  rownames(mat.phylo.sig) <- paste("Within",categories)
  return(as.data.frame(mat.phylo.sig))
}

rate.vs.tg <- list(
  kruskal_wallis_and_dunntest(formula = Net.diversification.rate..RPANDA. ~ Taxonomic.group,dat = data,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..RPANDA. ~ Taxonomic.group,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification.rate..RPANDA. ~ Taxonomic.group,dat = data,phy = tree),
  kruskal_wallis_and_dunntest(formula = Net.diversification.rate..CoMET. ~ Taxonomic.group,dat = data,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..CoMET. ~ Taxonomic.group,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification.rate..CoMET. ~ Taxonomic.group,dat = data,phy = tree),
  kruskal_wallis_and_dunntest(formula = Net.diversification.rate..ClaDS. ~ Taxonomic.group,dat = data,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..ClaDS. ~ Taxonomic.group,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification.rate..ClaDS. ~ Taxonomic.group,dat = data,phy = tree),
  kruskal_wallis_and_dunntest(formula = Net.diversification_Rate ~ Taxonomic.group,dat = data2,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_tip_rates(formula = Net.diversification.rate..ClaDS. ~ Taxonomic.group,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification_Rate ~ Taxonomic.group,dat = data2,phy = tree2)
  )

rate.vs.bo <- list(
  kruskal_wallis_and_dunntest(formula = Net.diversification.rate..RPANDA. ~ Biogeographic.origin,dat = data,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..RPANDA. ~ Biogeographic.origin,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification.rate..RPANDA. ~ Biogeographic.origin,dat = data,phy = tree),
  kruskal_wallis_and_dunntest(formula = Net.diversification.rate..CoMET. ~ Biogeographic.origin,dat = data,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..CoMET. ~ Biogeographic.origin,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification.rate..CoMET. ~ Biogeographic.origin,dat = data,phy = tree),
  kruskal_wallis_and_dunntest(formula = Net.diversification.rate..ClaDS. ~ Biogeographic.origin,dat = data,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..ClaDS. ~ Biogeographic.origin,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification.rate..ClaDS. ~ Biogeographic.origin,dat = data,phy = tree),
  kruskal_wallis_and_dunntest(formula = Net.diversification_Rate ~ Biogeographic.origin,dat = data2,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..ClaDS. ~ Biogeographic.origin,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification_Rate ~ Biogeographic.origin,dat = data2,phy = tree2)
)

rate.vs.bo <- list(
  kruskal_wallis_and_dunntest(formula = Net.diversification.rate..RPANDA. ~ Habitat,dat = data,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..RPANDA. ~ Habitat,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification.rate..RPANDA. ~ Habitat,dat = data,phy = tree),
  kruskal_wallis_and_dunntest(formula = Net.diversification.rate..CoMET. ~ Habitat,dat = data,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..CoMET. ~ Habitat,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification.rate..CoMET. ~ Habitat,dat = data,phy = tree),
  kruskal_wallis_and_dunntest(formula = Net.diversification.rate..ClaDS. ~ Habitat,dat = data,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..ClaDS. ~ Habitat,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification.rate..ClaDS. ~ Habitat,dat = data,phy = tree),
  kruskal_wallis_and_dunntest(formula = Net.diversification_Rate ~ Habitat,dat = data2,comp.name = "Diversification rate vs. Taxonomic group"),
  phyloANOVA_clade_rates(formula = Net.diversification.rate..ClaDS. ~ Habitat,phy = tree,df.pgls = data),
  phylosig_in_residuals(formula = Net.diversification_Rate ~ Habitat,dat = data2,phy = tree2)
)

rates_vs_clade_age <- list(phylosig_in_residuals(formula = Net.diversification.rate..RPANDA. ~ Crown.age..Mya.,dat = data,phy = tree),
                           phylosig_in_residuals(formula = Net.diversification.rate..CoMET. ~ Crown.age..Mya.,dat = data,phy = tree),
                           phylosig_in_residuals(formula = Net.diversification.rate..ClaDS. ~ Crown.age..Mya.,dat = data,phy = tree))


stats.rate.vs.tg <- list(anova_with_post_hoc(Net.diversification.rate..CoMET. ~  Taxonomic.group,data,comp.name = "Diversification rate vs. Taxonomic group"),
                      mk_pgls_plot(Net.diversification.rate..CoMET. ~  Taxonomic.group,T,tree,data,make.plot = FALSE),
                      phylosig_in_residuals(Net.diversification.rate..CoMET. ~  Taxonomic.group,data,tree),
                      phylosig_in_trait_within_categories_of_a_categorical_predictor(Net.diversification.rate..CoMET. ~  Taxonomic.group,data,tree))

stats.rate.vs.bo <- list(anova_with_post_hoc(Net.diversification.rate..CoMET. ~  Biogeographic.origin,data,comp.name = "Diversification rate vs. Biogeographic origin"),
                      mk_pgls_plot(Net.diversification.rate..CoMET. ~  Biogeographic.origin,T,tree,data,make.plot = FALSE),
                      phylosig_in_residuals(Net.diversification.rate..CoMET. ~  Biogeographic.origin,data,tree),
                      phylosig_in_trait_within_categories_of_a_categorical_predictor(Net.diversification.rate..CoMET. ~  Biogeographic.origin,data,tree))

stats.rate.vs.ht <- list(anova_with_post_hoc(Net.diversification.rate..CoMET. ~  Habitat,data,comp.name = "Diversification rate vs. Habitat type"),
                      mk_pgls_plot(Net.diversification.rate..CoMET. ~  Habitat,T,tree,data,make.plot = FALSE),
                      phylosig_in_residuals(Net.diversification.rate..CoMET. ~  Habitat,data,tree),
                      phylosig_in_trait_within_categories_of_a_categorical_predictor(Net.diversification.rate..CoMET. ~  Habitat,data,tree))

phyloanova.stats <- rbind("Taxonomic groups (PhyloANOVA)"=stats.rate.vs.tg[[2]]$phyloANOVA.results.phylANOVA,
                         "Biogeographic origins (PhyloANOVA)"=stats.rate.vs.bo[[2]]$phyloANOVA.results.phylANOVA,
                         "Habitat types (PhyloANOVA)"=stats.rate.vs.ht[[2]]$phyloANOVA.results.phylANOVA)

anova.stats <- rbind("Taxonomic groups (ANOVA)"=stats.rate.vs.tg[[1]]$ANOVA_stats,
                    "Biogeographic origins (ANOVA)"=stats.rate.vs.bo[[1]]$ANOVA_stats,
                    "Habitat types (ANOVA)"=stats.rate.vs.ht[[1]]$ANOVA_stats)

phyloanova.p.values <- c(stats.rate.vs.tg[[2]]$phyloANOVA.post.hoc,
                         stats.rate.vs.ht[[2]]$phyloANOVA.post.hoc)


anova.p.values <- c(stats.rate.vs.tg[[1]]$TukeyHSD_pairiwise_comparisons,
                    stats.rate.vs.ht[[1]]$TukeyHSD_pairiwise_comparisons)

phylogentic.signal <- rbind("Residuals\n(Rate vs. Taxonomic group)"=unlist(stats.rate.vs.tg[[3]][c("lambda","lambda.p.val","K","K.p.val")]),
                            stats.rate.vs.tg[[4]],
                            "Residuals\n(Rate vs. Biogeographic origin)"=unlist(stats.rate.vs.bo[[3]][c("lambda","lambda.p.val","K","K.p.val")]),
                            stats.rate.vs.bo[[4]],
                            "Residuals\n(Rate vs. Habitat type)"=unlist(stats.rate.vs.ht[[3]][c("lambda","lambda.p.val","K","K.p.val")]),
                            stats.rate.vs.ht[[4]]
                            )

#table.for.printing <- as.data.frame(p.vals.and.pagels.lambda)
table.for.printing1 <- as.data.frame(rbind(anova.stats,phyloanova.stats))
table.for.printing1 <- table.for.printing1[order(rownames(table.for.printing1),decreasing = T),]
table.for.printing1 <- table.for.printing1[c(2,1,6,5,4,3),]

table.for.printing2 <- as.data.frame(rbind(anova.p.values,phyloanova.p.values))
colnames(table.for.printing2) <- substr(colnames(table.for.printing2),27,nchar(colnames(table.for.printing2)))
colnames(table.for.printing2) <- str_replace(colnames(table.for.printing2),"-","-\n")
rownames(table.for.printing2) <- c("ANOVA-Tukey","PhyloANOVA-Tukey")

table.for.printing3 <- as.data.frame(phylogentic.signal)
inference.lambda <- inference.K <- vector(mode="character",length=nrow(table.for.printing3))

for(i in 1:nrow(table.for.printing3)){
  #i <- 3
  if(table.for.printing3$lambda.p.val[i]<0.05){
    inference.lambda[i] <- "significantly high\nphylogenetic signal"
  }else{
    if(round(table.for.printing3$lambda[i],1)<0.3){
      inference.lambda[i] <- "low phylogenetic signal\n(not significant)"
    }else{
      inference.lambda[i] <- "high phylogenetic signal\n(not significant)"
    }
  }
  
  if(table.for.printing3$K[i]<0.3){
    if(table.for.printing3$K.p.val[i]<0.05){
      inference.K[i] <- "significantly low\nphylogenetic signal"
    }else{
      inference.K[i] <- "low phylogenetic signal\n(not significant)"
    }
  }else{
    if(table.for.printing3$K.p.val[i]<0.05){
      inference.K[i] <- "significantly high\nphylogenetic signal"
    }else{
      inference.K[i] <- "high phylogenetic signal\n(not significant)"
    }
  }
}
table.for.printing3 <- cbind(table.for.printing3[,c(1,2)],inference.lambda,table.for.printing3[,c(3,4)],inference.K)
colnames(table.for.printing3) <- c("Pagel\'s \u03BB",
                                   "p-value for \u03BB",
                                   "Inference\n(from \u03BB)",
                                   "Blomberg\'s K",
                                   "p-value for K",
                                   "Inference\n(from K)")

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.45)),
  colhead = list(fg_params=list(cex = 0.5)),
  rowhead = list(fg_params=list(cex = 0.5)))

table.for.printing1 <- round(table.for.printing1,6)
table.for.printing2 <- round(table.for.printing2,6)
table.for.printing3 <- round(table.for.printing3,6)

stats.table <- tableGrob(table.for.printing1,theme = mytheme)
p.val.table <- tableGrob(table.for.printing2,theme = mytheme)
phylosig.table <- tableGrob(table.for.printing3,theme = mytheme)
valigned <- gtable_combine(stats.table,p.val.table,phylosig.table,along=2)
Cairo::CairoPDF("../graphs/tables_clade_rates.pdf",11.7,8.3,bg="transparent")
grid.arrange(valigned,nrow=1)
dev.off()

p1 <- ggplot(data,aes(x = Taxonomic.group, y = Net.diversification.rate..CoMET.,colour = Taxonomic.group))+
  xlab("Taxonomic group")+
  ylab("Net-diversification rate\n(events/lineage/Myr)")+
  geom_boxplot(color="black",outliers = F)+geom_jitter(alpha=0.4,size=3,stroke=0)+scale_color_manual(values = c("#00ff00","#000000","#0000ff"))+
  theme_classic()+theme(legend.position = "top",axis.title = element_text(size=10,face = "bold"),
                        axis.text = element_text(size=8))+guides(color=guide_legend(title="Taxonomic\ngroups"))

p2 <- ggplot(data,aes(x = Biogeographic.origin, y = Net.diversification.rate..CoMET.,colour = Biogeographic.origin))+
  xlab("Biogeographic origin")+
  ylab("Net-diversification rate\n(events/lineage/Myr)")+
  geom_boxplot(color="black",outliers = F)+geom_jitter(alpha=0.4,size=3,stroke=0)+scale_color_manual(values=c("grey10","grey40"))+
  theme_classic()+theme(legend.position = "top",axis.title = element_text(size=10,face = "bold"),
                        axis.text = element_text(size=8))+
  guides(color=guide_legend(title="Biogeographic\norigins"))

p3 <- ggplot(data,aes(x = Habitat, y = Net.diversification.rate..CoMET.,colour = Habitat))+
  xlab("Habitat")+
  ylab("Net-diversification rate\n(events/lineage/Myr)")+
  geom_boxplot(color="black",outliers = F)+geom_jitter(alpha=0.4,size=3,stroke=0)+scale_color_manual(values=c("grey10","grey40","grey70"))+
  theme_classic()+theme(legend.position = "top",axis.title = element_text(size=10,face = "bold"),
                        axis.text = element_text(size=8))+
  guides(color=guide_legend(title="Habitat\ntypes"))

p4 <- mk_pgls_plot(Net.diversification.rate..CoMET. ~  Crown.age..Mya.,is.predictor.categorical = F,phy=tree,df.pgls = data,xlab="Clade Age (Mya)",ylab="",
                   make.plot = TRUE,xlim.min = 0,ylim.min = 0)+theme_classic()+
  guides(color=guide_legend(title="Taxonomic\ngroups"))+theme(legend.position = "top",axis.title = element_text(size=10,face = "bold"),axis.text = element_text(size=8))


Cairo::CairoPDF("../graphs/clade_rate_vs_predictors.pdf",11.7,8.3,bg="transparent")
grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()



###################### TIP RATES ####################################

setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../data/"))

data2 <- read.csv("clads_tip_rates.csv")
data2 <- data2[-which(data2$Lineage=="Montecincla"),]
data2$Taxonomic.group <- data$Taxonomic.group[match(data2$Lineage,data$Lineage)]
data2$Biogeographic.origin <- data$Biogeographic.origin[match(data2$Lineage,data$Lineage)]
data2$Habitat <- data$Habitat[match(data2$Lineage,data$Lineage)]
data2$Crown.age..Mya. <- data$Crown.age..Mya.[match(data2$Lineage,data$Lineage)]
data2$No..of.species <- data$No..of.species[match(data2$Lineage,data$Lineage)]
tree2 <- read.tree("SUPER_TREE_with_all_tips_ultrametric.nwk")
tree2 <- drop.tip(tree2,tip = tree2$tip.label[grep(pattern = "Montecincla",tree2$tip.label)])

lineages <- data$Lineage

data2$tip_names <- rep(NA,nrow(data2))
for(i in 1:length(lineages)){
  #i <- 1
  ind.in.tree <- grep(lineages[i],tree2$tip.label)
  ind.in.matrix <- grep(lineages[i],data2$Lineage)
  data2$tip_names[ind.in.matrix] <- tree2$tip.label[ind.in.tree]
}



data2 <- data2[match(tree2$tip.label,data2$tip_names),] # order data in the same format as tip labels of tree
write.csv(data2,"clads_tip_rates_with_attributes.csv")


phyloANOVA_tip_rates <- function(formula,phy,df.pgls){
  # formula = Net.diversification_Rate ~ Taxonomic.group
  # phy = tree2
  # df.pgls <- data2
  
  Response <- as.character(formula)[2]
  Predictor <- as.character(formula)[3]
  
  colnames(df.pgls)[which(colnames(df.pgls) == Response)] <- "Response"
  colnames(df.pgls)[which(colnames(df.pgls) == Predictor)] <- "Predictor"
  
  
  df.pgls <- df.pgls[match(phy$tip.label,df.pgls$tip_names),] # order data in the same format as tip labels of tree

  tip.names <- df.pgls$tip_names
  rownames(df.pgls) <- tip.names
  
  Predictor <- as.factor(df.pgls$Predictor)
  names(Predictor) <- tip.names
  Response <- df.pgls$Response
  names(Response) <- tip.names
  # make variance covariance matrix from phylogeny
  cov.mat <- vcv.phylo(phy = phy)
  cat("\nPerforming GLS\n")
  phylo.anova <- lm.rrpp(Response~Predictor,
                         RRPP = T,
                         Cov = cov.mat,
                         print.progress = T,
                         verbose = T)
  
  sum.phylo.anova <- summary(phylo.anova)
  phylo.anova.table <- sum.phylo.anova$table[c("Df","SS","Residual SS","Rsq","F","Pr(>F)")]
  names(phylo.anova.table) <- c("df","SSR","SSE","MSR","F-value","p-value")
  
  
  if(length(unique(Predictor.mod))>2){
    cat("\nPerforming Post-Hoc\n")
    summary.pairwise.comp <- summary(pairwise(fit = phylo.anova,groups = Predictor,
                                              print.progress = T,verbose = T))
    
    phylANOVA.post.hoc <- convert.diag.matrix.to.vector(summary.pairwise.comp$pairwise.tables$P)
  }else{
    phylANOVA.post.hoc <- NA
  }
  cat("\nDone with Post-Hoc test\n")
  residuals <- residuals(phylo.anova)
  
  return(list("residuals.phyloANOVA"=residuals[,1],
              "phyloANOVA.results"=phylo.anova.table,
              "phyloANOVA.post.hoc"=phylANOVA.post.hoc))
  
}

# Phylogenetic signal in residuals
phylosig_in_residuals_tip_rates <- function(formula,dat,phy){
  # dat <- data2
  # phy <- tree2
  # formula <- Net.diversification_Rate ~ Biogeographic.origin
  lm.model <- lm(formula = formula,data=dat)
  residuals <- residuals(lm.model)
  names(residuals) <- dat$tip_names
  residuals <- residuals[match(phy$tip.label,names(residuals))]
  
  lambda.res <- phylosig(tree = phy,x = residuals,method = "lambda",test = T)
  K.res <- phylosig(tree = phy,x = residuals,method = "K",test = T,nsim = 1000)
  return(list("lambda"=lambda.res$lambda,"lambda.p.val"=lambda.res$P,"K"=K.res$K,"K.p.val"=K.res$P))
}


formula.tg <- Net.diversification_Rate ~ Taxonomic.group
lm.model.tg <- lm(formula = formula.tg,data=data2)
residuals.tg <- residuals(lm.model.tg)
names(residuals.tg) <- data2$tip_names
residuals.tg <- residuals.tg[match(tree2$tip.label,names(residuals.tg))]

formula.bo <- Net.diversification_Rate ~ Biogeographic.origin
lm.model.bo <- lm(formula = formula.bo,data=data2)
residuals.bo <- residuals(lm.model.bo)
names(residuals.bo) <- data2$tip_names
residuals.bo <- residuals.bo[match(tree2$tip.label,names(residuals.bo))]

formula.ht <- Net.diversification_Rate ~ Habitat
lm.model.ht <- lm(formula = formula.ht,data=data2)
residuals.ht <- residuals(lm.model.ht)
names(residuals.ht) <- data2$tip_names
residuals.ht <- residuals.ht[match(tree2$tip.label,names(residuals.ht))]



df <- as.data.frame(cbind(data2$Net.diversification_Rate,residuals.tg,residuals.bo,residuals.ht))
head(df)
colnames(df) <- c("Rates","Residuals for TG","Residuals for BO","Residuals for HT")
data3 <- data2[match(rownames(df),data2$tip_names),]




# Phylogenetic signal in traits within groups
phylosig_in_trait_within_categories_of_a_categorical_predictor_tip_rates <- function(formula,dat,phy){
  # dat <- data2;phy <- tree2
  # formula <- Net.diversification_Rate ~ Taxonomic.group
  Response <- as.character(formula)[2]
  Predictor <- as.character(formula)[3]
  
  categories <- unique(dat[,Predictor])
  mat.phylo.sig <- matrix(nrow=length(categories),ncol=4)
  rownames(mat.phylo.sig) <- categories
  colnames(mat.phylo.sig) <- c("lambda","lambda.p.val","K","K.p.val")
  for(g in categories){
    # g <- categories[1]
    dat.g <- dat[which(dat[,Predictor] == g),]
    phy.g <- keep.tip(phy,dat.g$tip_names)
    dat.g <- dat.g[match(phy.g$tip.label,dat.g$tip_names),]
    trait.g <- dat.g[,Response]
    names(trait.g) <- dat.g$tip_names
    mat.phylo.sig[g,c("lambda","lambda.p.val")] <- unlist(phylosig(tree = phy.g, x = trait.g,method = "lambda",test = T)[c("lambda","P")])
    mat.phylo.sig[g,c("K","K.p.val")] <- unlist(phylosig(tree = phy.g, x = trait.g,method = "K",test = T,nsim = 1000)[c("K","P")])
  }
  rownames(mat.phylo.sig) <- paste("Within",categories)
  return(as.data.frame(mat.phylo.sig))
}

mk_phylotrait_maps <- function(phy,df){
  p4d <- phylo4d(tree2,df)
  
  cols <- c(rep(c("blue","red"),16),"blue")
  
  names(cols) <- unique(data2$Lineage)
  
  cols <- cols[match(data2$Lineage,names(cols))]
  cols <- cols[match(data3$Lineage,names(cols))]
  names(cols) <- data3$tip_names
  pdf("../graphs/tree_with_residuals_tip_rates.pdf",11.7,8.3,bg="transparent")
  # phylosignal::barplot.phylo4d(p4d)
  # phylosignal::gridplot(p4d, tree.type = "fan", tree.ratio = 0.5,
  #                       show.trait = FALSE, show.tip = T,
  #                       cell.col = terrain.colors(100))
  
  phylosignal::barplot.phylo4d(p4d, show.tip = F, tree.open.angle = 30,
                               bar.lwd=1,bar.col = cols,
                               tree.ratio = 0.6,trait.cex = 0.8,tip.cex=0.1)
  dev.off()
}


stats.rate.vs.tg.tip_rates <- list(anova_with_post_hoc(Net.diversification_Rate ~  Taxonomic.group,data2,comp.name = "Diversification rate vs. Taxonomic group"),
                         phyloANOVA_tip_rates(Net.diversification_Rate ~  Taxonomic.group,tree2,data2),
                         phylosig_in_residuals_tip_rates(Net.diversification_Rate ~  Taxonomic.group,data2,tree2),
                         phylosig_in_trait_within_categories_of_a_categorical_predictor_tip_rates(Net.diversification_Rate ~  Taxonomic.group,data2,tree2))

stats.rate.vs.bo.tip_rates <- list(anova_with_post_hoc(Net.diversification_Rate ~  Biogeographic.origin,data2,comp.name = "Diversification rate vs. Biogeographic origin"),
                                   phyloANOVA_tip_rates(Net.diversification_Rate ~  Biogeographic.origin,tree2,data2),
                                   phylosig_in_residuals_tip_rates(Net.diversification_Rate ~  Biogeographic.origin,data2,tree2),
                                   phylosig_in_trait_within_categories_of_a_categorical_predictor_tip_rates(Net.diversification_Rate ~  Biogeographic.origin,data2,tree2))

stats.rate.vs.ht.tip_rates <- list(anova_with_post_hoc(Net.diversification_Rate ~  Habitat,data2,comp.name = "Diversification rate vs. Habitat"),
                                   phyloANOVA_tip_rates(Net.diversification_Rate ~  Habitat,tree2,data2),
                                   phylosig_in_residuals_tip_rates(Net.diversification_Rate ~  Habitat,data2,tree2),
                                   phylosig_in_trait_within_categories_of_a_categorical_predictor_tip_rates(Net.diversification_Rate ~  Habitat,data2,tree2))

phyloanova.stats.tip_rates <- rbind("Taxonomic groups (PhyloANOVA)"=stats.rate.vs.tg.tip_rates[[2]]$phyloANOVA.results,
                          "Biogeographic origins (PhyloANOVA)"=stats.rate.vs.bo.tip_rates[[2]]$phyloANOVA.results,
                          "Habitat types (PhyloANOVA)"=stats.rate.vs.ht.tip_rates[[2]]$phyloANOVA.results)

anova.stats.tip_rates <- rbind("Taxonomic groups (ANOVA)"=stats.rate.vs.tg.tip_rates[[1]]$ANOVA_stats,
                     "Biogeographic origins (ANOVA)"=stats.rate.vs.bo.tip_rates[[1]]$ANOVA_stats,
                     "Habitat types (ANOVA)"=stats.rate.vs.ht.tip_rates[[1]]$ANOVA_stats)

phyloanova.p.values.tip_rates <- c(stats.rate.vs.tg.tip_rates[[2]]$phyloANOVA.post.hoc,
                         stats.rate.vs.ht.tip_rates[[2]]$phyloANOVA.post.hoc)


anova.p.values.tip_rates <- c(stats.rate.vs.tg.tip_rates[[1]]$TukeyHSD_pairiwise_comparisons,
                    stats.rate.vs.ht.tip_rates[[1]]$TukeyHSD_pairiwise_comparisons)

phylogentic.signal.tip_rates <- rbind("Residuals\n(Rate vs. Taxonomic group)"=unlist(stats.rate.vs.tg.tip_rates[[3]][c("lambda","lambda.p.val","K","K.p.val")]),
                            stats.rate.vs.tg.tip_rates[[4]],
                            "Residuals\n(Rate vs. Biogeographic origin)"=unlist(stats.rate.vs.bo.tip_rates[[3]][c("lambda","lambda.p.val","K","K.p.val")]),
                            stats.rate.vs.bo.tip_rates[[4]],
                            "Residuals\n(Rate vs. Habitat type)"=unlist(stats.rate.vs.ht.tip_rates[[3]][c("lambda","lambda.p.val","K","K.p.val")]),
                            stats.rate.vs.ht.tip_rates[[4]])

#table.for.printing <- as.data.frame(p.vals.and.pagels.lambda)
table.for.printing1.tip_rates <- as.data.frame(rbind(anova.stats.tip_rates,phyloanova.stats.tip_rates))
table.for.printing1.tip_rates <- table.for.printing1.tip_rates[order(rownames(table.for.printing1.tip_rates),decreasing = T),]
table.for.printing1.tip_rates <- table.for.printing1.tip_rates[c(2,1,6,5,4,3),]

table.for.printing2.tip_rates <- as.data.frame(rbind(anova.p.values.tip_rates,phyloanova.p.values.tip_rates))
colnames(table.for.printing2.tip_rates) <- substr(colnames(table.for.printing2.tip_rates),27,nchar(colnames(table.for.printing2.tip_rates)))
colnames(table.for.printing2.tip_rates) <- str_replace(colnames(table.for.printing2.tip_rates),"-","-\n")
rownames(table.for.printing2.tip_rates) <- c("ANOVA-Tukey (p-value)","PhyloANOVA-Tukey (p-value)")

table.for.printing3.tip_rates <- as.data.frame(phylogentic.signal.tip_rates)
inference.lambda <- inference.K <- vector(mode="character",length=nrow(table.for.printing3.tip_rates))

for(i in 1:nrow(table.for.printing3.tip_rates)){
  #i <- 1
  if(table.for.printing3.tip_rates$lambda.p.val[i]<0.05){
    inference.lambda[i] <- "significantly high\nphylogenetic signal"
  }else{
    if(round(table.for.printing3.tip_rates$lambda[i],1)<0.3){
      inference.lambda[i] <- "low phylogenetic signal\n(not significant)"
    }else{
      inference.lambda[i] <- "high phylogenetic signal\n(not significant)"
    }
  }
  
  if(table.for.printing3.tip_rates$K[i]<0.3){
    if(table.for.printing3.tip_rates$K.p.val[i]<0.05){
      inference.K[i] <- "significantly low\nphylogenetic signal"
    }else{
      inference.K[i] <- "low phylogenetic signal\n(not significant)"
    }
  }else{
    if(table.for.printing3.tip_rates$K.p.val[i]<0.05){
      inference.K[i] <- "significantly high\nphylogenetic signal"
    }else{
      inference.K[i] <- "high phylogenetic signal\n(not significant)"
    }
  }
}
table.for.printing3.tip_rates <- cbind(table.for.printing3.tip_rates[,c(1,2)],
                                       inference.lambda,
                                       table.for.printing3.tip_rates[,c(3,4)],
                                       inference.K)
colnames(table.for.printing3.tip_rates) <- c("Pagel\'s \u03BB",
                                   "p-value for \u03BB",
                                   "Inference\n(from \u03BB)",
                                   "Blomberg\'s K",
                                   "p-value for K",
                                   "Inference\n(from K)")

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.5)),
  colhead = list(fg_params=list(cex = 0.6,fontface=c(rep("bold",3),rep("bold.italic",3)))),
  rowhead = list(fg_params=list(cex = 0.6)))
mytheme2 <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.5)),
  colhead = list(fg_params=list(cex = 0.6)),
  rowhead = list(fg_params=list(cex = 0.6)))

table.for.printing1.tip_rates <- round(table.for.printing1.tip_rates,6)
table.for.printing2.tip_rates <- round(table.for.printing2.tip_rates,6)
table.for.printing3.tip_rates <- round(table.for.printing3.tip_rates,6)
table.for.printing1.tip_rates <- table.for.printing1.tip_rates[,-1]
# paste(table.for.printing3.tip_rates$`Pagel's λ`[c(1,5,8)],
#        paste0("(",table.for.printing3.tip_rates$`p-value for λ`[c(1,5,8)],")"),sep = "\n")

values.to.plug <- paste(round(table.for.printing3.tip_rates$`Pagel's λ`[c(1,5,8)],6),
                                                                 paste0("(<0.001*)"),sep = "\n")

table.for.printing1.tip_rates[,"Pagel's \u03BB\n(p-value for \u03BB)"] <- rep(NA,nrow(table.for.printing1.tip_rates))
table.for.printing1.tip_rates[c(2,4,6),"Pagel's \u03BB\n(p-value for \u03BB)"] <- values.to.plug

table.for.printing1.tip_rates$`p-value`[which(table.for.printing1.tip_rates$`p-value`<0.001)] <- "<0.001*"

stats.table.tip_rates <- tableGrob(table.for.printing1.tip_rates,theme = mytheme)
stats.table.tip_rates <- gtable::gtable_add_grob(stats.table.tip_rates,
                                         grobs=grid::rectGrob(gp=grid::gpar(fill=NA,lwd=2)),
                                         t=1,b=nrow(stats.table.tip_rates),l=1,
                                         r=ncol(stats.table.tip_rates))
# stats.table.tip_rates <- gtable::gtable_add_grob(stats.table.tip_rates,
#                                                  grobs=grid::rectGrob(gp=grid::gpar(fill=NA,lwd=2)),
#                                                  t=1,l=1,
#                                                  r=ncol(stats.table.tip_rates))

p.val.table.tip_rates <- tableGrob(table.for.printing2.tip_rates,theme = mytheme2)
p.val.table.tip_rates <- gtable::gtable_add_grob(p.val.table.tip_rates,
                                                 grobs=grid::rectGrob(gp=grid::gpar(fill=NA,lwd=2)),
                                                 t=1,b=nrow(p.val.table.tip_rates),l=1,
                                                 r=ncol(p.val.table.tip_rates))
# p.val.table.tip_rates <- gtable::gtable_add_grob(p.val.table.tip_rates,
#                                                  grobs=grid::rectGrob(gp=grid::gpar(fill=NA,lwd=2)),
#                                                  t=1,l=1,
#                                                  r=ncol(p.val.table.tip_rates))
#phylosig.table.tip_rates <- tableGrob(table.for.printing3.tip_rates,theme = mytheme)
valigned <- gtable_combine(stats.table.tip_rates,p.val.table.tip_rates,along=2)
Cairo::CairoPDF("../graphs/tables_tip_rates.pdf",8.3,11.7,bg="transparent")
grid.arrange(valigned,nrow=2)
dev.off()

p1.tip_rate <- ggplot(data2,aes(x = Taxonomic.group, y = Net.diversification_Rate,fill = Taxonomic.group,color=Taxonomic.group))+
  xlab("Taxonomic group")+
  ylab("Net-diversification rate\n(events/lineage/Myr)")+
  # geom_violin(aes(fill = Taxonomic.group, fill = after_scale(colorspace::lighten(fill, .5))),
  #             size = 0.2, bw = .05) +
  geom_boxplot(width=0.9, alpha=0.3,outliers = F,lwd=0.4) +
  stat_dots(side = "both",color=NA,alpha=0.3,dotsize=1.07)+
  scale_fill_manual(values = c("#77dd77","#898989","#00c7e7"))+
  scale_color_manual(values = c("#77dd77","#898989","#00c7e7"))+
  theme_classic()+
  theme(legend.position = "none",axis.title = element_text(size=20,face = "bold"),
                        axis.text = element_text(size=15))

p2.tip_rate <- ggplot(data2,aes(x = Biogeographic.origin, y = Net.diversification_Rate,fill = Biogeographic.origin,color = Biogeographic.origin))+
  xlab("Biogeographic origin")+
  ylab("")+
  # geom_violin(aes(fill = Biogeographic.origin, fill = after_scale(colorspace::lighten(fill, .5))),
  #             size = 0.2, bw = .05) +
  geom_boxplot(width=0.9, alpha=0.3,outliers = F,lwd=0.4) +
  stat_dots(side = "both",color=NA,alpha=0.3,dotsize=1.07)+
  scale_fill_manual(values = c("grey10","grey40"))+
  scale_color_manual(values = c("grey10","grey40"))+
  theme_classic()+
  theme(legend.position = "none",axis.title = element_text(size=20,face = "bold"),
        axis.text = element_text(size=15))
  

p3.tip_rate <- ggplot(data2,aes(x = Habitat, y = Net.diversification_Rate,fill = Habitat, color = Habitat))+
  xlab("Habitat types")+
  ylab("Net-diversification rate\n(events/lineage/Myr)")+
  # geom_violin(aes(fill = Habitat, fill = after_scale(colorspace::lighten(fill, .5))),
  #             size = 0.2, bw = .05) +
  geom_boxplot(width=0.9, alpha=0.3,outliers = F,lwd=0.4) +
  stat_dots(side = "both",color=NA,alpha=0.3,dotsize=1.07)+
  scale_fill_manual(values = c("grey10","grey40","grey70"))+
  scale_color_manual(values = c("grey10","grey40","grey70"))+
  theme_classic()+
  theme(legend.position = "none",axis.title = element_text(size=20,face = "bold"),
        axis.text = element_text(size=15))

data$log_clade_age  <- log(data$Crown.age..Mya.)
p4 <- mk_pgls_plot(Net.diversification.rate..CoMET. ~  Crown.age..Mya.,is.predictor.categorical = F,
                   phy=tree,df.pgls = data,xlab="Clade Age (Mya)",
                   ylab="",
                   make.plot = TRUE,xlim.min=2,ylim.min=0)+
  theme(legend.position = c(0.585,0.85),legend.direction = "vertical",
        legend.box = "horizontal",legend.background = element_blank(),legend.text = element_text(size = 10),
        legend.box.background = element_rect(colour = "black",size = 0.1),
        panel.border = element_blank(),axis.title = element_text(size=20,face = "bold"),
        axis.text = element_text(size=15))


Cairo::CairoPDF("../graphs/tip_and_clade_rate_vs_predictors.pdf",11.7,8.3,bg="transparent")
grid.arrange(p1.tip_rate,p2.tip_rate,p3.tip_rate,p4,nrow=2)
dev.off()
  
