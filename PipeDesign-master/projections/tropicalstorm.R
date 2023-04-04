# Bayesian Inference via MCMC
#Initialize Data
rm(list=ls())
rm(list=ls())
options(scipen=9999)
main_path=getwd()

infile<-paste(main_path,"/dataset/AMS_ellicott.txt",sep='')
big<-scan(infile,skip=1,list(date=0,depth=0)) 
depth<-big$depth 
datset<-depth*25.4/24

filename<-paste(main_path,"/dataset/Tropicalstorms.txt",sep='')
big<-scan(filename,skip=101,list(date=0,storms=0,hurricanes=0,majorhurricanes=0,ace=0))
year<-big$date
TS_hist<-big$storms
time_hist<-year
# re-normalize
ind_norm <- which(time_hist==1951):which(time_hist==2018) 
tempset <- (TS_hist - mean(TS_hist[ind_norm]))/sd(TS_hist[ind_norm])

#####projections
infile<-paste(main_path,"/dataset/AMS_mdrtempGFDL.txt",sep='')
big<-scan(infile,skip=0,list(date=0,depth=0))
temperature_proj<-big$depth
time_proj<-big$date
year_proj<-time_proj
ibeg<-which(year_proj==1951)
iend<-which(year_proj==2018)
temp_proj=temperature_proj-mean(temperature_proj[ibeg:iend])
ind_proj <- which(year_proj==1951):which(year_proj==2018)
temp_proj=(temp_proj - mean(temp_proj[ind_proj]))/sd(temp_proj[ind_proj])
temp_proj=temp_proj[which(year_proj==2020):which(year_proj==2090)]

# (1) Stationary
source(paste(main_path,"/SourceCode/Prior2SourceStat.R",sep=""))
# Initial Conditions
start<-c(rep(1,3)) # Start Value
errvect<-c(0.0625, 0.03125, 0.06255) # Proposal Distribution SD
uvect<-c(rep(0,3)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(100,3))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors
#MCMC
run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2")
# Diagnostics
# Run Diagnostics Again
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=6000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
save(run2,file=paste(main_path,"/TSRESULT/stat_widenorm_run2_TS.RData",sep=""))
cred.table
save(cred.table,file=paste(main_path,"/TSRESULT/stat_widenorm_CI_TS.RData",sep=""))
mu=run2$finmat[,1]
sigma=exp(run2$finmat[,2])
xi=run2$finmat[,3]
mu_chain <-mu[(length(xi)-40000+1):length(xi)]
sigma_chain <- sigma[(length(sigma)-40000+1):length(sigma)]
xi_chain <- xi[(length(xi)-40000+1):length(xi)]
retint<-benreturn(100,mu_chain,sigma_chain,xi_chain)
save(mu_chain,xi_chain,sigma_chain,file=paste(main_path,"/TSRESULT/stat_widenorm_param_TS.RData",sep=""))
save(retint,file=paste(main_path,"/TSRESULT/stat_widenorm_rtnlevel.RData",sep=""))

# (2) Mu Non-Stationary
source(paste(main_path,"/SourceCode/Prior2SourceMu.R",sep=""))
# Initial Conditions
start<-c(rep(1,4)) # Start Value
errvect<-c(0.0625 ,0.03125, 0.0625 ,0.0625) # SD for Proposal 
uvect<-c(rep(0,4)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(100,4))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors
#MCMC
run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Diagnostics
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=10000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
save(run2,file=paste(main_path,"/TSRESULT/Mu_nonstat_widenorm_run2_RCM.RData",sep=""))
cred.table
save(cred.table,file=paste(main_path,"/TSRESULT/Mu_nonstat_widenorm_CI_RCM.RData",sep=""))
mu0=run2$finmat[,1]
amu=run2$finmat[,4]
mu=mu0*(1+amu*tempset[68])
sigma=exp(run2$finmat[,2])
xi=run2$finmat[,3]
mu_chain <-mu[(length(xi)-40000+1):length(xi)]
mu0_chain <-mu0[(length(xi)-40000+1):length(xi)]
amu0_chain <-amu[(length(xi)-40000+1):length(xi)]
slope<- mu0_chain*amu0_chain
sigma_chain <- sigma[(length(sigma)-40000+1):length(sigma)]
xi_chain <- xi[(length(xi)-40000+1):length(xi)]
retint<-benreturn(100,mu_chain,sigma_chain,xi_chain)
save(mu_chain,xi_chain,sigma_chain,file=paste(main_path,"/TSRESULT/Mu_nonstat_widenorm_param_TS.RData",sep=""))
save(retint,file=paste(main_path,"/TSRESULT/Mu_nonstat_widenorm_rtnlevel.RData",sep=""))
save(mu0_chain,amu0_chain,slope,file=paste(main_path,"/TSRESULT/slope_nonstat_widenorm_param_TS.RData",sep=""))
