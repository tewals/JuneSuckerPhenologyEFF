#######################################################
# June Sucker Migration Timing Analysis
#   (Walsworth et al. 2024, Ecology of Freshwater Fish)
#
#   Code to run migration timing model in JAGS
#    and to produce basic summary figures
#
#   - Includes code for both the normal run timing
#      model and the gamma run timing model
#     - Manuscript focuses on normal run timing
#       model due to its superior fit
#
# Contact: timothy.walsworth@usu.edu
#----------------------------------------------------

####Setup Libraries and import data####
library(jagsUI) 
library(postpack)
library(R2jags)
library(runjags)
library(rjags)


# import the csv with migration data 
dat<-read.csv("JS.minmax.csv")

#### Defining data values nad vectors to feed into JAGS model

# Vector of total hits per day = Nt (all days and years)
Nt<- dat$Ntotal
# Vector of year = year (all days and years)
year<- (dat$Year-2007) #dirty but works year 1 = 2008, 2=2009 etc

# Vector of DOY = doy (all days and years)
doy.g<- dat$DOY_shift  # used in gamma model
doy.n<- dat$DOY # used in normal model

# Make a data frame of just the years with directional data 
dat.us<-dat[dat$Year>2017,]
dat.us[is.na(dat.us)]<-0

#' Vector of integers linking year indices between 
#' short and full datasets = ylink
#  - A vector of years with up/dn data = yearsup
ylink<- (dat.us$Year - 2007) 

# Vector of doys with up/dn data = doyup
doyup.g<- dat.us$DOY_shift
doyup.n<- dat.us$DOY

# Integer total number of observations across all years = nobsall
nobsall<- length(dat$Ntotal) 

# Integer number of obs with up/dn data = nobsup
nobsup<- length(dat.us$Ntotal)

#Add for the binomProb
Nup.obs<- dat.us$Nupstream #number of upstream observations 
Ndn.obs<-dat.us$Ndownstream

# Integer number of years up/dn data = nyrup
nyrup<- length(unique(yearsup))

####Model Setup and Run#### 
mod.data<- list(Nt=Nt, 
                year=year, 
                doy=doy.n, 
                Nup.obs=Nup.obs, 
                Ndn.obs=Ndn.obs,
                doyup=doyup.n, 
                nobsall=nobsall,
                nobsup=nobsup,
                nyr=nyr, 
                ylink=ylink)

# Normal Run Timing Model 
# - R working directory needs to be set to folder containing
#   this .bug file
mod.loc<-"JS_MigrationTiming_Model_Normal_JAGScode.bug"  

# #MCMC settings 
nc<- 3 #number of chains

# Specify Model Parameters for JAGS
mod.params<- c('peak.date', 'lsig.run', 'lrun.yr', 
               'mu.peak', 'sig.peak', 'lmu.run', 'lsig.run.hp',
               'wt', 'thetat','llN','llNup','llNdn',
               'mu.therealsig.run','sig.therealsig.run',
               'd','shape.d','rate.d') 

# run the model in parallel
start.time<-Sys.time()
mod.out<-jags.parallel(data=mod.data,parameters.to.save = mod.params,n.chains=nc,n.iter=300000,n.burnin=120000,
              n.thin = 500,model.file = mod.loc)
Sys.time()-start.time

plot(mod.out) # Check model estimates
#traceplot(mod.out)

# Extract MCMC posteriors
mcout<-mod.out$BUGSoutput$sims.list

# Extract Posterior distributions for model parameters
post.peak.date<- (as.matrix(mcout$peak.date)) # Peak date
post.run.yr<- exp(as.matrix(mcout$lrun.yr)) # Run size
post.sig.run<- exp(as.matrix(mcout$lsig.run)) # Among-individual variation (sd)
post.d<-as.matrix(mcout$d) # Residence time


####
# Predict upstream and downstream migrations per day for all years
#  given parameter estimates from model

dys<-seq(1,365) # Days of year to simulate run size across for each year

pred.dbd<-pred.rbd<-array(0,dim=c(nyr, length(dys) ,nrow(post.peak.date))) # Storage for predicted runs
d10n<-d90n<-matrix(0,nrow=nrow(post.peak.date),ncol=nyr)
for(i in 1:nrow(post.peak.date)){
  for(j in 1:nyr){
    
    # Predict Upstream migrants per day
    pred.rbd[j,,i]<-post.run.yr[i,j]/(post.sig.run[i,j]*sqrt(2*pi))*
      exp(-((dys-post.peak.date[i,j])^2/(2*post.sig.run[i,j]^2)))
    
    # Predict downstream migrants per day
    pred.dbd[j,,i]<-post.run.yr[i,j]/(post.sig.run[i,j]*sqrt(2*pi))*
      exp(-(((dys-post.d[i,j])-post.peak.date[i,j])^2/(2*post.sig.run[i,j]^2)))
    
    # Calculate correction factor psi
    psi<-sum(1/(post.sig.run[i,j]*sqrt(2*pi))*
               exp(-((dys-post.peak.date[i,j])^2/(2*post.sig.run[i,j]^2))))
    
    # Apply correction factor to Upstream and downstream movements
    pred.rbd[j,,i]<-pred.rbd[j,,i]/psi
    pred.dbd[j,,i]<-pred.dbd[j,,i]/psi
    
    # Calculate the day when 10% of run has entered the river
    d10n[i,j]<-qnorm(.10,mean=post.peak.date[i,j],sd=post.sig.run[i,j])
    
    # Calculate the day when 90% of the run has left the river
    d90n[i,j]<-qnorm(.90,mean=post.peak.date[i,j]+post.d[i,j],sd=post.sig.run[i,j])
    
  }
  
  
}

pred.tbd<-pred.rbd+pred.dbd # total movements by day

# Create storage matrices for predicted quantiles
pred75u<-pred25u<-predhiu<-predlou<-predmedu<-matrix(0,nrow=length(dys),ncol=nyr)
pred75d<-pred25d<-predhid<-predlod<-predmedd<-predmedu
pred75t<-pred25t<-predhit<-predlot<-predmedt<-predmedu

# Calculate Quantiles of predicted migrants per day across all years
# - rows are day of year, columns are year
for(i in 1:nyr){
  for(j in 1:length(dys)){
    # Upstream movements
    predmedu[j,i]<-median(pred.rbd[i,j,])
    predhiu[j,i]<-quantile(pred.rbd[i,j,],probs=0.975)
    predlou[j,i]<-quantile(pred.rbd[i,j,],probs=0.025)
    pred75u[j,i]<-quantile(pred.rbd[i,j,],probs=0.75)
    pred25u[j,i]<-quantile(pred.rbd[i,j,],probs=0.25)
    
    # Downstream movements
    predmedd[j,i]<-median(pred.dbd[i,j,])
    predhid[j,i]<-quantile(pred.dbd[i,j,],probs=0.975)
    predlod[j,i]<-quantile(pred.dbd[i,j,],probs=0.025)
    pred75d[j,i]<-quantile(pred.dbd[i,j,],probs=0.75)
    pred25d[j,i]<-quantile(pred.dbd[i,j,],probs=0.25)
    
    # Total movements
    predmedt[j,i]<-median(pred.tbd[i,j,])
    predhit[j,i]<-quantile(pred.tbd[i,j,],probs=0.975)
    predlot[j,i]<-quantile(pred.tbd[i,j,],probs=0.025)
    pred75t[j,i]<-quantile(pred.tbd[i,j,],probs=0.75)
    pred25t[j,i]<-quantile(pred.tbd[i,j,],probs=0.25)
  }
}

####
# Plot model predictions with observed data in a series of 
#  2x2 panel plots

layout(matrix(c(1,2,
                3,4),nrow=2,ncol=2,byrow=T))
par(mar=c(3,4.5,1,1))
for(i in 2008:2022){
  
  plot(dat$DOY[dat$Year == i],dat$Ntotal[dat$Year == i],pch=16,col="darkgreen",xlab="",
       ylab="Number",xlim=c(70,220),ylim=c(0,1.05*max(c(predhit[,i-2007],dat$Ntotal[dat$Year==i]))))
  lines(predmedt[,i-2007],lwd=2,col="darkgreen")
  if(i%in%c(2018:2022)){
    points(dat.us$DOY[dat.us$Year == i],dat.us$Nupstream[dat.us$Year == i],pch=16,col="dodgerblue")
    points(dat.us$DOY[dat.us$Year == i],dat.us$Ndownstream[dat.us$Year == i],pch=16,col="darkorange")
    
  }
  lines(predmedu[,i-2007],lwd=2,col="dodgerblue")
  
  lines(predmedd[,i-2007],lwd=2,col="darkorange")
  legend("topright",bty="n",legend=i)

}

# Read in data to convert Julian Dates for plot labels
doyconv<-read.csv("DOYConversion.csv",header=T)
months<-c("January","February","March","April","May","June","July",
          "August","September","October","November","December")

# Plot posterior distributions for peak date, among individual variation
#  and residence time across years
layout(matrix(c(1,2,3),nrow=1,ncol=3))
yoff<-seq(14,0,by=-1)
mday<-rep(0,nyr)

# For peak date
plot(NA,NA,ylim=c(0,nyr),xlim=c(100,170),xlab="",
     ylab="",yaxt="n",bty="n",yaxs="i",xaxt="n")
axis(side=2,at=seq(0,14),labels=seq(2022,2008,by=-1),las=1)
axis(side=1,at=doyconv$DOY[doyconv$Day==1],labels=paste(months,"1"),las=2)
axis(side=1,at=doyconv$DOY[doyconv$Day==15],labels=paste(months,"15"),las=2)
abline(h=seq(0,14),col="grey")
for(i in 2008:2022){
  ds<-density(mcout$peak.date[,i-2007])
  mday[i-2007]<-(ds$x[which.max(ds$y)])
  lines(ds$x,0.95*(ds$y)/max(ds$y)+yoff[i-2007],lwd=2,
        col=ifelse(i<2018,"grey","black"))

}

# For Among individual variation
plot(NA,NA,ylim=c(0,nyr),xlim=c(0,50),xlab="Days",
     ylab="",yaxt="n",bty="n",yaxs="i",xaxs="i")
axis(side=2,at=seq(0,14),labels=seq(2022,2008,by=-1),las=1)

abline(h=seq(0,14),col="grey")
for(i in 2008:2022){
  ds<-density(exp(mcout$lsig.run[,i-2007]))
  lines(ds$x,0.95*(ds$y)/max(ds$y)+yoff[i-2007],lwd=2,
        col=ifelse(i<2018,"grey","black"))

}

# For residence time
plot(NA,NA,ylim=c(0,nyr),xlim=c(0,40),xlab="Days",
     ylab="",yaxt="n",bty="n",yaxs="i",xaxs="i")
axis(side=2,at=seq(0,14),labels=seq(2022,2008,by=-1),las=1)

abline(h=seq(0,14),col="grey")
for(i in 2008:2022){
  
  ds<-density(mcout$d[,i-2007])
  lines(ds$x,0.95*(ds$y)/max(ds$y)+yoff[i-2007],lwd=2,
        col=ifelse(i<2018,"grey","black"))
  
}
mtext(outer=T,side=3,line=-3.7,at=c(.18,.5,.85),text=c("Peak Date","SD of\nUpstream Migration","Residence Time"))


###############
### Gamma run timing model
###############


####Model Setup and Run#### 
mod.data<- list(Nt=Nt, 
                year=year, 
                doy=doy.g, 
                Nup.obs=Nup.obs, 
                Ndn.obs=Ndn.obs,
                doyup=doyup.g, 
                nobsall=nobsall,
                nobsup=nobsup,
                nyr=nyr, 
                ylink=ylink,
                mindoy.yr=mindoy.yr)


# #MCMC settings 
nc<- 3 #number of chains

mod.loc<-"JSMigration_Analysis_updown_Gamma_stdPriors.bug"  # Gamma Run Timing Model 

#parameters for updown normal
mod.params<- c('peak.date', 'lsig.run','lrun.yr','mu.peak', 'sig.peak',
               'lmu.run', 'lsig.run.hp','wt', 'thetat',
               'llN','llNup','llNdn','mu.therealsig.run',
               'sig.therealsig.run',
               'd','shape.d','rate.d','Ndn.pred2') 

# Run model in parallel
start.time<-Sys.time()
mod.out<-jags.parallel(data=mod.data,parameters.to.save = mod.params,
                       n.chains=nc,n.iter=300000,n.burnin=120000,
                       n.thin = 500,model.file = mod.loc)
Sys.time()-start.time

plot(mod.out) # check model outputs
#traceplot(mod.out)

# Extract MCMC posteriors
mcout<-mod.out$BUGSoutput$sims.list


# Extract Posterior distributions for model parameters
post.peak.date<- exp(as.matrix(mcout$peak.date)) # peak date
post.run.yr<- exp(as.matrix(mcout$lrun.yr)) # Run size
post.sig.run<- exp(as.matrix(mcout$lsig.run)) # Among-individual variation
post.d<-as.matrix(mcout$d) # residence time


####
# Predict upstream and downstream migrations per day for all years
#  given parameter estimates from model

dys<-seq(1,365) # Days of year to simulate run size across for each year

# Create storage objects for predicted runs
pred.dbd<-pred.rbd<-array(0,dim=c(nyr, length(dys) ,nrow(post.peak.date)))
d10g<-d90g<-matrix(0,nrow=nrow(post.peak.date),ncol=nyr)

for(i in 1:nrow(post.peak.date)){ # For each posterior sample
  for(j in 1:nyr){ # for each year
    
    # Sample from posteriors
    peakuse<-post.peak.date[i,j]
    peakused<-post.peak.date[i,j]+post.d[i,j]
    
    # Calculate gamma parameters from peak, among-individual variation parameters
    rate<-(peakuse+sqrt(peakuse^2+4*post.sig.run[i,j]^2))/(2*post.sig.run[i,j]^2)
    shape<-peakuse*rate+1
    
    ratedn<-(peakused+sqrt(peakused^2+4*post.sig.run[i,j]^2))/(2*post.sig.run[i,j]^2)
    shapedn<-peakused*ratedn+1
    
    # Predict upstream movements per day
    pred.rbd[j,,i]<-((post.run.yr[i,j])*(rate^shape)/exp(lgamma(shape))*dys^(shape-1)*
                        exp(-rate*dys))
    
    # Predict downstream movements per day
    pred.dbd[j,,i]<-((post.run.yr[i,j])*(ratedn^shapedn)/exp(lgamma(shapedn))*dys^(shapedn-1)*
                       exp(-ratedn*dys))
  
    # Calculate correction factor psi for upstream moves
    psi<-sum((rate^shape)/exp(lgamma(shape))*dys^(shape-1)*
               exp(-rate*dys))
    
    # Calculate correction factor psi for downstream moves
    psid<-sum((ratedn^shapedn)/exp(lgamma(shapedn))*dys^(shapedn-1)*
               exp(-ratedn*dys))
    
    # Correct for psi
    pred.rbd[j,,i]<-pred.rbd[j,,i]/psi
    pred.dbd[j,,i]<-pred.dbd[j,,i]/psid
    
    # Determine day in which 10% of run has entered river
    d10g[i,j]<-qnorm(.10,mean=post.peak.date[i,j],sd=post.sig.run[i,j])
    
    # Calculate day in which 90 % of run has left river
    d90g[i,j]<-qnorm(.90,mean=post.peak.date[i,j]+post.d[i,j],sd=post.sig.run[i,j])
    
  }
  
  
}

# Calculate total movements by day
pred.tbd<-pred.rbd+pred.dbd

# Create storage for quantiles of daily predicted movements
pred75u<-pred25u<-predhiu<-predlou<-predmedu<-matrix(0,nrow=length(dys),ncol=nyr)
pred75d<-pred25d<-predhid<-predlod<-predmedd<-predmedu
pred75t<-pred25t<-predhit<-predlot<-predmedt<-predmedu

# Calculate quantiles of daily predicted movements across all years
for(i in 1:nyr){
  for(j in 1:length(dys)){
    # Upstream
    predmedu[j,i]<-median(pred.rbd[i,j,],na.rm=T)
    predhiu[j,i]<-quantile(pred.rbd[i,j,],probs=0.975,na.rm=T)
    predlou[j,i]<-quantile(pred.rbd[i,j,],probs=0.025,na.rm=T)
    pred75u[j,i]<-quantile(pred.rbd[i,j,],probs=0.75,na.rm=T)
    pred25u[j,i]<-quantile(pred.rbd[i,j,],probs=0.25,na.rm=T)
    
    # Downstream
    predmedd[j,i]<-median(pred.dbd[i,j,],na.rm=T)
    predhid[j,i]<-quantile(pred.dbd[i,j,],probs=0.975,na.rm=T)
    predlod[j,i]<-quantile(pred.dbd[i,j,],probs=0.025,na.rm=T)
    pred75d[j,i]<-quantile(pred.dbd[i,j,],probs=0.75,na.rm=T)
    pred25d[j,i]<-quantile(pred.dbd[i,j,],probs=0.25,na.rm=T)
    
    # Total
    predmedt[j,i]<-median(pred.tbd[i,j,],na.rm=T)
    predhit[j,i]<-quantile(pred.tbd[i,j,],probs=0.975,na.rm=T)
    predlot[j,i]<-quantile(pred.tbd[i,j,],probs=0.025,na.rm=T)
    pred75t[j,i]<-quantile(pred.tbd[i,j,],probs=0.75,na.rm=T)
    pred25t[j,i]<-quantile(pred.tbd[i,j,],probs=0.25,na.rm=T)
  }
}

# Plot median model predictions against observed data for each year
#  in a series of 2x2 panel figures
layout(matrix(c(1,2,
                3,4),nrow=2,ncol=2,byrow=T))
par(mar=c(3,4.5,1,1))

for(i in 2008:2022){

  plot(dat$DOY[dat$Year == i],dat$Ntotal[dat$Year == i],pch=16,col="darkgreen",xlab="",
       ylab="Number",xlim=c(80,200),
       ylim=c(0,1.05*max(c(predhit[,i-2007],dat$Ntotal[dat$Year==i]))))
 
  lines(dys+mindoy.yr[i-2007]-1, predmedt[,i-2007],lwd=2,col="darkgreen")
  
  if(i%in%c(2018:2022)){
    points(dat.us$DOY[dat.us$Year == i],
           dat.us$Nupstream[dat.us$Year == i],pch=16,col="dodgerblue")
    points(dat.us$DOY[dat.us$Year == i],
           dat.us$Ndownstream[dat.us$Year == i],pch=16,col="darkorange")
    
  }
  
  lines(dys+mindoy.yr[i-2007]-1,predmedu[,i-2007],lwd=2,col="dodgerblue")
  lines(dys+mindoy.yr[i-2007]-1,predmedd[,i-2007],lwd=2,col="darkorange")
  legend("topright",bty="n",legend=i)

}

# Read in data to convert Julian Dates for plot labels
doyconv<-read.csv("DOYConversion.csv",header=T)
months<-c("January","February","March","April","May","June","July",
          "August","September","October","November","December")

# Plot posterior distributions of peak date, among individual variation
#   and residence time for each year, stacked to allow comparison.
layout(matrix(c(1,2,3),nrow=1,ncol=3))
yoff<-seq(14,0,by=-1)

# For peak date
plot(NA,NA,ylim=c(0,nyr),xlim=c(90,136),xlab="",
     ylab="",yaxt="n",bty="n",yaxs="i",xaxt="n")
axis(side=2,at=seq(0,14),labels=seq(2022,2008,by=-1),las=1)
axis(side=1,at=doyconv$DOY[doyconv$Day==1],labels=paste(months,"1"),las=2)
axis(side=1,at=doyconv$DOY[doyconv$Day==15],labels=paste(months,"15"),las=2)
abline(h=seq(0,14),col="grey")
for(i in 2008:2022){
  ds<-density(mcout$peak.date[,i-2007])
  lines(ds$x+mindoy.yr[i-2007],0.95*(ds$y)/max(ds$y)+yoff[i-2007],lwd=2,
        col=ifelse(i<2018,"black","dodgerblue"))
}

# For among-individual variation
plot(NA,NA,ylim=c(0,nyr),xlim=c(0,50),xlab="Days",
     ylab="",yaxt="n",bty="n",yaxs="i",xaxs="i")
axis(side=2,at=seq(0,14),labels=seq(2022,2008,by=-1),las=1)
abline(h=seq(0,14),col="grey")
for(i in 2008:2022){
  ds<-density(exp(mcout$lsig.run[,i-2007]))
  lines(ds$x,0.95*(ds$y)/max(ds$y)+yoff[i-2007],lwd=2,
        col=ifelse(i<2018,"black","dodgerblue"))
}

# For residence time
plot(NA,NA,ylim=c(0,nyr),xlim=c(0,100),xlab="Days",
     ylab="",yaxt="n",bty="n",yaxs="i",xaxs="i")
axis(side=2,at=seq(0,14),labels=seq(2022,2008,by=-1),las=1)
abline(h=seq(0,14),col="grey")
for(i in 2008:2022){
  ds<-density(d90[,i-2007]-d10[i-2007])
  lines(ds$x,0.95*(ds$y)/max(ds$y)+yoff[i-2007],lwd=2,
        col=ifelse(i<2018,"black","dodgerblue"))
}

mtext(outer=T,side=3,line=-3.7,at=c(.18,.5,.85),text=c("Peak Date","SD of\nUpstream Migration","Residence Time"))
