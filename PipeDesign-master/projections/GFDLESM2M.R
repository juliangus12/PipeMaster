#Initialize Data
rm(list=ls())
options(scipen=9999)
main_path=getwd()
infile<-paste(main_path,"/dataset/maca/AMS_GFDL_ESM2M.txt",sep='')
big<-scan(infile,skip=64,list(date=0,depth=0)) 
depthrcp<-big$depth/24 
datset<-depthrcp[7:77] 
##Projection
library(ncdf4)
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
tempset=temp_proj[which(year_proj==2020):which(year_proj==2090)]
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
save(run2,file=paste(main_path,"/GFDLESM2M/Mu_nonstat_widenorm_run2_RCM.RData",sep=""))
cred.table
save(cred.table,file=paste(main_path,"/GFDLESM2M/Mu_nonstat_widenorm_CI_RCM.RData",sep=""))
mu0=run2$finmat[,1]
amu=run2$finmat[,4]
sigma=exp(run2$finmat[,2])
xi=run2$finmat[,3]
g=0
gg=2019
 for (i in 1:71){
 g=g+1
 gg=gg+1
 muproj=mu0*(1+amu*tempset[i])
 muproj_chain <-muproj[(length(xi)-40000+1):length(xi)]
 sigma_chain <- sigma[(length(sigma)-40000+1):length(sigma)]
 xi_chain <- xi[(length(xi)-40000+1):length(xi)]
save(muproj_chain,xi_chain,sigma_chain,file=paste(main_path,'/GFDLESM2M/nonstat_param',gg,'.RData',sep=""))}
