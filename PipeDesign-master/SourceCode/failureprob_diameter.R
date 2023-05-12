
library(Hmisc)
library(gplots)
library(plotrix)
library(ggplot2)
library(evir)
library(plyr)

rm(list=ls())
options(scipen=9999)
main_path=getwd()


## function to estimate the return level from GEV distribution
myreturnlevel <- function(t,mu,sigma,xi){
  x<-qgev(p=1-(1/t),xi=xi,mu=mu,sigma=sigma)
  return(x)}
#--------------------------------------------------------------

####Historical_AMS_dataset
infile<-paste(main_path,"/dataset/AMS_ellicott.txt",sep='')
big<-scan(infile,skip=1,list(date=0,depth=0))
depth<-big$depth 
datset<-depth*25.4/24

####NOAA atlas 14
infile<-paste(main_path,"/dataset/atlas14_Ellicott.txt",sep='')
big<-scan(infile,skip=5,list(rtn=0,noaa=0,noaalow=0,noaaup=0))
rtn<-big$rtn
noaa<-big$noaa
noaalow<-big$noaalow
noaaup<-big$noaaup

######stat
paramstat<-load(paste(main_path,"/MDRRESULT/stat_widenorm_param_MDR.RData",sep=""))
paramtime1<-mu_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
st100<-myreturnlevel(100, paramtime1, (paramtime2), paramtime3)

####Test for different pipe diameters
dia1= 0.3048  #12inches
dia2=0.381    #15 in
dia3= 0.4572  #18 in
dia4=0.5334   #21 inches  ####NOAA DIAMETER
dia5=0.6096   #24 inches
dia6=0.6858   #27 inches
dia7=0.762    #30 inches
dia8=0.9144   #36 inches
dia9=1.0668  #42 inches
dia10=1.2192  #48 inches
#####assume different runoff coefficients 
runoffcoeff<-c(0.45,0.55,0.65,0.75,0.85)

####for each available pipe diameter, compute failure probaility under different rainfall projections
###nonstat
g<-0
for (i in 2020:2090){
infile=paste(main_path,'/MDRRESULT/nonstat_param',i,'.RData',sep="")
load(infile)
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
nst100<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
stat100<-c(runoffcoeff[1]*st100,runoffcoeff[2]*st100,runoffcoeff[3]*st100,runoffcoeff[4]*st100,runoffcoeff[5]*st100)  
nstat100<-c(runoffcoeff[1]*nst100,runoffcoeff[2]*nst100,runoffcoeff[3]*nst100,runoffcoeff[4]*nst100,runoffcoeff[5]*nst100)
nNOAA100<-c(runoffcoeff[3]*9)
area=0.25
####compute flow load capacity#####
Qstat100<-0.278*stat100*area  
Qnstat100<-0.278*nstat100*area 
Qnoaa100<-0.278*nNOAA100*area  
slope<-0.01   #assumption
mann<-0.013   #assumption
####design diameter using noaa estimates
designdia<-((Qnoaa100*mann)/(0.31*slope^0.5))^(3/8)   
designdia
###compute pipe's flow capacity
Qcap<-c((0.31/mann)*dia4^(8/3)*slope^(1/2))  
####compute failure probability#####
GQstat100<-Qcap-Qstat100          
Gnstat100<-Qcap-Qnstat100
Hstat100<-count(GQstat100<0)
Hnstat100<-count(Gnstat100<0)
length<-length(GQstat100)
Pstat100<-Hstat100[2,2]/length   
g<-g+1
Pnstat100[g]<-Hnstat100[2,2]/length 
write.table(Pnstat100,file=paste("MDRRESULT_dia4.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}
