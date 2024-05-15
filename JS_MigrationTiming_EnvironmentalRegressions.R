#----------------------------------------------------
# June Sucker Migration Timing Analysis
#   (Walsworth et al. 2024, Ecology of Freshwater Fish)
#
#   Environmental Predictors of Migration Phenology
#    Regressions
#
# Contact: timothy.walsworth@usu.edu
#----------------------------------------------------

#' Load Migration Timing Model Posterior MCMC samples
load("JS_MigrationTiming_Model_Posteriors.rdata")

#' Use the postpack package to extract the posterior
#'  distribution for peak migration date (peak.date),
#'  among individual variation (lsig.run), and
#'  residence time (d)
library(postpack)
post.peak.date<- (as.matrix(mod.out$BUGSoutput$sims.list$peak.date))
post.sig.run <- (as.matrix(exp(mod.out$BUGSoutput$sims.list$lsig.run)))
post.res.time <- (as.matrix(mod.out$BUGSoutput$sims.list$d))


#' Read in environmental data for each year
env.dat <- read.csv("JS_MigrationTiming_EnvironmentalData.csv",header=T)

#' Create vectors for PC1, PC2, and PC3 to incorporate
#'  in linear regressions
pc1 <- env.dat$PC1 # flow volume
pc2 <- env.dat$PC2 # flow timing
pc3 <- env.dat$PC3 # Temp


#' Function to run linear regressions for a given response variable
lrFUN <- function(post.data, resp.id){
  #' Create an array to hold model output from our 4 alternate models
  lm.out <- array(0, dim=c(nrow(post.data),12,10)) 
  colnames(lm.out) <- c("intercept","PC1","PC2","PC3",
                        "PC1:PC2","PC1:PC3","PC2:PC3","PC!:PC2:PC3",
                        "rsq","AIC","resp","pval")
  
  #' Models to run
  #' 1: response ~ PC1
  #' 2: response ~ PC2
  #' 3: response ~ PC1 + PC2
  #' 4: response ~ PC1*PC2 == PC1 + PC2 + PC1:PC2
  #' 5: response ~ PC3
  #' 6: response ~ PC1 + PC3
  #' 7: response ~ PC2 + PC3
  #' 8: response ~ PC1 + PC2 + PC3
  
  
  #' Run a for loop to run linear regressions for each
  #'  draw from the posterior (first dimension of the array)
  #'  for each of the alternate models 
  #'  (third dimension of output array).
  #' Each model has 15 data points for the 15 years of data
  
  for(i in 1:nrow(lm.out)){ # For each posterior sample
    
    # Fit linear regression model with only PC1 as predictor
    mod1 <- lm(post.data[i,] ~ pc1) 

    lm.out[i,1,1] <- coef(mod1)[1] # Store intercept estimate
    lm.out[i,2,1] <- coef(mod1)[2] # Store PC1 slope estimate
    lm.out[i,9,1] <- summary(mod1)$r.squared # Store r2
    lm.out[i,10,1] <- AIC(mod1) # Store AIC
    lm.out[i,11,1] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod1)$fstatistic # Extract Fstatistics
    lm.out[i,12,1] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    # Fit linear regression model with only PC2 as predictor
    mod2 <- lm(post.data[i,] ~ pc2)
    
    lm.out[i,1,2] <- coef(mod2)[1] # Store intercept estimate
    lm.out[i,3,2] <- coef(mod2)[2] # Store PC2 slope estimate
    lm.out[i,9,2] <- summary(mod2)$r.squared # Store r2
    lm.out[i,10,2] <- AIC(mod2) # Store AIC
    lm.out[i,11,2] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod2)$fstatistic # Extract Fstatistics
    lm.out[i,12,2] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    # Fit linear regression model with only PC3 as predictor
    mod3 <- lm(post.data[i,] ~ pc3)
    
    lm.out[i,1,3] <- coef(mod3)[1] # Store intercept estimate
    lm.out[i,4,3] <- coef(mod3)[2] # Store PC3 slope estimate
    lm.out[i,9,3] <- summary(mod3)$r.squared # Store r2
    lm.out[i,10,3] <- AIC(mod3) # Store AIC
    lm.out[i,11,3] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod3)$fstatistic # Extract Fstatistics
    lm.out[i,12,3] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    # Fit linear regression model with additive PC1 and PC2 as predictors
    mod4 <- lm(post.data[i,] ~ pc1 + pc2)
    
    lm.out[i,1,4] <- coef(mod4)[1] # Store intercept estimate
    lm.out[i,2,4] <- coef(mod4)[2] # Store PC1 slope estimate
    lm.out[i,3,4] <- coef(mod4)[3] # Store PC2 slope estimate
    lm.out[i,9,4] <- summary(mod4)$r.squared # Store r2
    lm.out[i,10,4] <- AIC(mod4) # Store AIC
    lm.out[i,11,4] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod4)$fstatistic # Extract Fstatistics
    lm.out[i,12,4] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    # Fit linear regression model with additive PC1 and PC3 as predictors
    mod5 <- lm(post.data[i,] ~ pc1+pc3)
    
    lm.out[i,1,5] <- coef(mod5)[1] # Store intercept estimate
    lm.out[i,2,5] <- coef(mod5)[2] # Store PC1 slope estimate
    lm.out[i,4,5] <- coef(mod5)[3] # Store PC3 slope estimate
    lm.out[i,9,5] <- summary(mod5)$r.squared # Store r2
    lm.out[i,10,5] <- AIC(mod5) # Store AIC
    lm.out[i,11,5] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod5)$fstatistic # Extract Fstatistics
    lm.out[i,12,5] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    # Fit linear regression model with additive PC2 and PC3 as predictors
    mod6 <- lm(post.data[i,] ~ pc2+pc3)
    
    lm.out[i,1,6] <- coef(mod6)[1] # Store intercept estimate
    lm.out[i,3,6] <- coef(mod6)[2] # Store PC2 slope estimate
    lm.out[i,4,6] <- coef(mod6)[3] # Store PC3 slope estimate
    lm.out[i,9,6] <- summary(mod6)$r.squared # Store r2
    lm.out[i,10,6] <- AIC(mod6) # Store AIC
    lm.out[i,11,6] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod6)$fstatistic # Extract Fstatistics
    lm.out[i,12,6] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    #' Fit linear regression model with additive PC1, PC2,
    #'  and PC3 as predictors
    mod7 <- lm(post.data[i,] ~ pc1+pc2+pc3)
    
    lm.out[i,1,7] <- coef(mod7)[1] # Store intercept estimate
    lm.out[i,2,7] <- coef(mod7)[2] # Store PC1 slope estimate
    lm.out[i,3,7] <- coef(mod7)[3] # Store PC2 slope estimate
    lm.out[i,4,7] <- coef(mod7)[4] # Store PC3 slope estimate
    lm.out[i,9,7] <- summary(mod7)$r.squared # Store r2
    lm.out[i,10,7] <- AIC(mod7) # Store AIC
    lm.out[i,11,7] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod7)$fstatistic # Extract Fstatistics
    lm.out[i,12,7] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    #' Fit linear regression model with interaction 
    #' between PC1 and PC2
    mod8 <- lm(post.data[i,] ~ pc1*pc2)
    
    lm.out[i,1,8] <- coef(mod8)[1] # Store intercept estimate
    lm.out[i,2,8] <- coef(mod8)[2] # Store PC1 slope estimate
    lm.out[i,3,8] <- coef(mod8)[3] # Store PC2 slope estimate
    lm.out[i,5,8] <- coef(mod8)[4] # Store interaction estimate
    lm.out[i,9,8] <- summary(mod8)$r.squared # Store r2
    lm.out[i,10,8] <- AIC(mod8) # Store AIC
    lm.out[i,11,8] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod8)$fstatistic # Extract Fstatistics
    lm.out[i,12,8] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    #' Fit linear regression model with interaction 
    #' between PC1 and PC3
    mod9 <- lm(post.data[i,] ~ pc1*pc3)
    
    lm.out[i,1,9] <- coef(mod9)[1] # Store intercept estimate
    lm.out[i,2,9] <- coef(mod9)[2] # Store PC1 slope estimate
    lm.out[i,4,9] <- coef(mod9)[3] # Store PC3 slope estimate
    lm.out[i,6,9] <- coef(mod9)[4] # Store interaction estimate
    lm.out[i,9,9] <- summary(mod9)$r.squared # Store r2
    lm.out[i,10,9] <- AIC(mod9) # Store AIC
    lm.out[i,11,9] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod9)$fstatistic # Extract Fstatistics
    lm.out[i,12,9] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    #' Fit linear regression model with interaction 
    #' between PC2 and PC3
    mod10 <- lm(post.data[i,] ~ pc2*pc3)
    
    lm.out[i,1,10] <- coef(mod10)[1] # Store intercept estimate
    lm.out[i,3,10] <- coef(mod10)[2] # Store PC2 slope estimate
    lm.out[i,4,10] <- coef(mod10)[3] # Store PC3 slope estimate
    lm.out[i,7,10] <- coef(mod10)[4] # Store interaction estimate
    lm.out[i,9,10] <- summary(mod10)$r.squared # Store r2
    lm.out[i,10,10] <- AIC(mod10) # Store AIC
    lm.out[i,11,10] <- resp.id # Store specified response ID (not used)
    fst<-summary(mod10)$fstatistic # Extract Fstatistics
    lm.out[i,12,10] <-pf(fst[1],fst[2],fst[3],lower.tail=FALSE) # Store p-value
    
    
  }
  
  #' Calculate mean R2 value across regressions for
  #'   all posterior draws
  
  mean.R2 <- rep(NA,10)
  for(i in 1:10){
    mean.R2[i] <- mean(lm.out[,9,i])
  }
  
  names(mean.R2) <- c("mod1","mod2","mod3","mod4","mod5","mod6","mod7",
                      "mod8","mod9","mod10")
  mean.R2 <- round(mean.R2,3)
  
  
  #' Calculate mean AIC value across regressions for
  #'   all posterior draws
  mean.AIC <- rep(NA,10)
  
  for(i in 1:10){
    mean.AIC[i] <- mean(lm.out[,10,i])
  }
  
  #' Calculate delta AIC from means for each alternative model
  deltaAIC <- mean.AIC-min(mean.AIC)
  names(deltaAIC) <- c("mod1","mod2","mod3","mod4","mod5","mod6","mod7",
                       "mod8","mod9","mod10")
  
  #' Return a list with all model outputs, mean R2, and delta AIC
  #'   values
  return(list(lm.out,mean.R2,deltaAIC))
}

#' Run linear regression models for peak migration date estimates
pDOYpeak<-lrFUN(post.data=post.peak.date,resp.id = 1)

#' Extract summary data from model with lowest deltaAIC value

#' Quantiles of R2 estimates
DOYpeak.r2.q<-quantile(pDOYpeak[[1]][,9,which.min(pDOYpeak[[3]])],probs=c(.025,.05,.25,.5,.75,.95,.975))

#' Quantiles of p-values
DOYpeak.pv.q<-quantile(pDOYpeak[[1]][,12,which.min(pDOYpeak[[3]])],probs=c(.025,.05,.25,.5,.75,.95,.975))

#' Proportion of p-values less than 0.05
DOYpeak.pv.less05<-sum(pDOYpeak[[1]][,12,which.min(pDOYpeak[[3]])]<0.05)/length(pDOYpeak[[1]][,12,which.min(pDOYpeak[[3]])])

#' Determine which Coefficients are included in the best 
#' fitting model for later extraction by identifying column
#' numbers with parameter estimates
sigpreds<-which((colSums(abs(pDOYpeak[[1]][,c(2:8),which.min(pDOYpeak[[3]])])))>0)+1


#' Run linear regression models for among-individual
#'  variation estimates
pDOYsd<-lrFUN(post.data=post.sig.run,resp.id = 1)

#' Extract summary data from model with lowest deltaAIC value

#' Quantiles of R2 estimates
DOYsd.r2.q<-quantile(pDOYsd[[1]][,9,which.min(pDOYsd[[3]])],probs=c(.025,.05,.25,.5,.75,.95,.975))

#' Quantiles of p-values
DOYsd.pv.q<-quantile(pDOYsd[[1]][,12,which.min(pDOYsd[[3]])],probs=c(.025,.05,.25,.5,.75,.95,.975))

#' Proportion of p-values less than 0.05
DOYsd.pv.less05<-sum(pDOYsd[[1]][,12,which.min(pDOYsd[[3]])]<0.05)/length(pDOYsd[[1]][,12,which.min(pDOYsd[[3]])])

#' Determine which Coefficients are included in the best 
#' fitting model for later extraction by identifying column
#' numbers with parameter estimates
sigpreds.sd<-which((colSums(abs(pDOYsd[[1]][,c(2:8),which.min(pDOYsd[[3]])])))>0)+1


#' Run linear regression models for residence time estimates
pDOYd<-lrFUN(post.data=post.res.time,resp.id = 1)

#' Extract summary data from model with lowest deltaAIC value

#' Quantiles of R2 estimates
DOYd.r2.q<-quantile(pDOYd[[1]][,9,which.min(pDOYd[[3]])],probs=c(.025,.05,.25,.5,.75,.95,.975))

#' Quantiles of p-values
DOYd.pv.q<-quantile(pDOYd[[1]][,12,which.min(pDOYd[[3]])],probs=c(.025,.05,.25,.5,.75,.95,.975))

#' Proportion of p-values less than 0.05
DOYd.pv.less05<-sum(pDOYd[[1]][,12,which.min(pDOYd[[3]])]<0.05)/length(pDOYd[[1]][,12,which.min(pDOYd[[3]])])

#' Determine which Coefficients are included in the best 
#' fitting model for later extraction by identifying column
#' numbers with parameter estimates
sigpreds.d<-which((colSums(abs(pDOYd[[1]][,c(2:8),which.min(pDOYd[[3]])])))>0)+1


#' Create list of regression parameter estimate outputs
#' for best fitting models for each migration timing
#' parameter
reg.out<-list(pDOYpeak[[1]][,c(1,sigpreds),which.min(pDOYpeak[[3]])],
              pDOYsd[[1]][,c(1,sigpreds.sd),which.min(pDOYsd[[3]])],
              pDOYd[[1]][,c(1,sigpreds.d),which.min(pDOYd[[3]])])

#' Create a list of summary stat quantiles for each migration
#' timing parameter
regressionstatsout<-list(reg.out,DOYpeak.r2.q,DOYpeak.pv.q,DOYpeak.pv.less05,
                         DOYsd.r2.q,DOYsd.pv.q,DOYsd.pv.less05,
                         DOYd.r2.q,DOYd.pv.q,DOYd.pv.less05)


