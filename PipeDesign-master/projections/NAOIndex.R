# Bayesian Inference via MCMC
rm(list=ls())
options(scipen=9999)
main_path=getwd()

infile<-paste(main_path,"/dataset/AMS_ellicott.txt",sep='')
big<-scan(infile,skip=1,list(date=0,depth=0))
depth<-big$depth 
datset<-depth*25.4/24

####nao_data
nao_dat <- read.table(paste(main_path,"/dataset/NAO.txt",sep=''))
colnames(nao_dat) <- c('year','jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec','ann')
ibeg <- which(nao_dat['year']==1951)
iend <- max(which(nao_dat[,'ann']!=-99.99))
time_hist <- nao_dat[ibeg:iend, 'year']
# get DJF means
nao_hist <- rep(-999, length(time_hist))
for (y in 1:length(time_hist)) {
  nao_hist[y] <- mean( c(nao_dat$jan[ibeg+y-1], nao_dat$feb[ibeg+y-1], nao_dat$mar[ibeg+y-1]) )
}
time_hist <- nao_dat[ibeg:iend, 'year']
# re-normalize
ind_norm <- which(time_hist==1951):which(time_hist==2018)
tempset <- (nao_hist - mean(nao_hist[ind_norm]))/sd(nao_hist[ind_norm])
######Projection
library(ncdf4)
ncdata<-nc_open(paste(main_path,"/dataset/DMIEH5_SRA1B_4_MM_psl.1-1200.nc",sep=''))
ncvar_get(ncdata)
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
burnin=10000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
cred.table
save(run2,file=paste(main_path,"/NAOIndex/stat_widenorm_run2_NAO.RData",sep=""))
save(cred.table,file=paste(main_path,"/NAOIndex/stat_widenorm_CI_NAO.RData",sep=""))
mu=run2$finmat[,1]
sigma=exp(run2$finmat[,2])
xi=run2$finmat[,3]
mu_chain <-mu[(length(xi)-40000+1):length(xi)]
sigma_chain <- sigma[(length(sigma)-40000+1):length(sigma)]
xi_chain <- xi[(length(xi)-40000+1):length(xi)]
retint<-benreturn(100,mu_chain,sigma_chain,xi_chain)
save(mu_chain,xi_chain,sigma_chain,file=paste(main_path,"/NAOIndex/stat_widenorm_param_NAO.RData",sep=""))
save(retint,file=paste(main_path,"/NAOIndex/stat_widenorm_rtnlevel.RData",sep=""))
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
burnin=6000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
save(run2,file=paste(main_path,"/NAOIndex/Mu_nonstat_widenorm_run2_NAO.RData",sep=""))
cred.table
save(cred.table,file=paste(main_path,"/NAOIndex/Mu_nonstat_widenorm_CI_NAO.RData",sep=""))
mu0=run2$finmat[,1]
amu=run2$finmat[,4]
mu=mu0*(1+amu*tempset[68])
sigma=exp(run2$finmat[,2])
xi=run2$finmat[,3]
mu_chain <-mu[(length(xi)-40000+1):length(xi)]
sigma_chain <- sigma[(length(sigma)-40000+1):length(sigma)]
xi_chain <- xi[(length(xi)-40000+1):length(xi)]
retint<-benreturn(100,mu_chain,sigma_chain,xi_chain)
save(mu_chain,xi_chain,sigma_chain,file=paste(main_path,"/NAOIndex/Mu_nonstat_widenorm_param_NAO.RData",sep=""))
save(retint,file=paste(main_path,"/NAOIndex/Mu_nonstat_widenorm_rtnlevel.RData",sep=""))
