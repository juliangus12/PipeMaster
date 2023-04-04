# Load libraries and data requidarkred to run this code

library(Hmisc)
library(gplots)
library(plotrix)
library(ggplot2)
library(evir)
main_path=getwd()
## function to estimate the return level from GEV distribution
myreturnlevel <- function(t,mu,sigma,xi){
  x<-qgev(p=1-(1/t),xi=xi,mu=mu,sigma=sigma)
  return(x)}
#--------------------------------------------------------------
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

infile<-paste(main_path,"/dataset/AMS_ellicott.txt",sep='')
big<-scan(infile,skip=1,list(date=0,depth=0))
depth<-big$depth 
datset<-depth*25.4/24

####NOAA atlas 14
infile<-paste(main_path,"/dataset/atlas14_Ellicott.txt",sep='')
big<-scan(infile,skip=1,list(rtn=0,noaa=0,noaalow=0,noaaup=0))
rtn<-big$rtn
noaa<-big$noaa
noaalow<-big$noaalow
noaaup<-big$noaaup

######stationary
paramstat<-load(paste(main_path,"/GTRESULT/stat_widenorm_param_GT.RData",sep=''))
paramtime1<-mu_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
MCrlstat <- sapply(1:length(paramtime1), function(x) {
  myreturnlevel((1.1:100.1), paramtime1[x], (paramtime2[x]), paramtime3[x])
})
MCrlstatmean <- apply(MCrlstat, 1, mean)
lowerstat <- sapply(1.1:100.1, function (x)
{  quantile(MCrlstat[x,], 0.05)  })
upperstat <- sapply(1.1:100.1, function (x)
{  quantile(MCrlstat[x,], 0.95) })

####Globaltemp
paramGT<-load(paste(main_path,"/GTRESULT/Mu_nonstat_widenorm_param_GT.RData",sep=''))
paramtime1<-mu_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
MCrlGT <- sapply(1:length(paramtime1), function(x) {
  myreturnlevel((1.1:100.1), paramtime1[x], (paramtime2[x]), paramtime3[x])
})
MCrlGTmean <- apply(MCrlGT, 1, mean)
lowerGT <- sapply(1.1:100.1, function (x)
{  quantile(MCrlGT[x,], 0.05)  })
upperGT <- sapply(1.1:100.1, function (x)
{  quantile(MCrlGT[x,], 0.95) })

####Localtemp
paramLT<-load(paste(main_path,"/LTRESULT/Mu_nonstat_widenorm_param_LT.RData",sep=''))
paramtime1<-mu_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
MCrlLT <- sapply(1:length(paramtime1), function(x) {
  myreturnlevel((1.1:100.1), paramtime1[x], (paramtime2[x]), paramtime3[x])
})
MCrlLTmean <- apply(MCrlLT, 1, mean)
lowerLT <- sapply(1.1:100.1, function (x)
{  quantile(MCrlLT[x,], 0.05)  })
upperLT <- sapply(1.1:100.1, function (x)
{  quantile(MCrlLT[x,], 0.95) })

#####MDRtemp
paramMDR<-load(paste(main_path,"/MDRRESULT/Mu_nonstat_widenorm_param_MDR.RData",sep=''))
paramtime1<-mu_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
MCrlMDR <- sapply(1:length(paramtime1), function(x) {
  myreturnlevel((1.1:100.1), paramtime1[x], (paramtime2[x]), paramtime3[x])
})
MCrlMDRmean <- apply(MCrlMDR, 1, mean)
lowerMDR <- sapply(1.1:100.1, function (x)
{  quantile(MCrlMDR[x,], 0.05)  })
upperMDR <- sapply(1.1:100.1, function (x)
{  quantile(MCrlMDR[x,], 0.95) })

######NAO Index
paramNAO<-load(paste(main_path,"/NAOIndex/Mu_nonstat_widenorm_param_NAO.RData",sep=''))
paramtime1<-mu_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
MCrlNAO <- sapply(1:length(paramtime1), function(x) {
  myreturnlevel((1.1:100.1), paramtime1[x], (paramtime2[x]), paramtime3[x])
})
MCrlNAOmean <- apply(MCrlNAO, 1, mean)
lowerNAO <- sapply(1.1:100.1, function (x)
{  quantile(MCrlNAO[x,], 0.05)  })
upperNAO <- sapply(1.1:100.1, function (x)
{  quantile(MCrlNAO[x,], 0.95) })

#####Tropicalstorm
paramTS<-load(paste(main_path,"/TSRESULT/Mu_nonstat_widenorm_param_TS.RData",sep=''))
paramtime1<-mu_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
MCrlTS <- sapply(1:length(paramtime1), function(x) {
  myreturnlevel((1.1:100.1), paramtime1[x], (paramtime2[x]), paramtime3[x])
})
MCrlTSmean <- apply(MCrlTS, 1, mean)
lowerTS <- sapply(1.1:100.1, function (x)
{  quantile(MCrlTS[x,], 0.05)  })
upperTS <- sapply(1.1:100.1, function (x)
{  quantile(MCrlTS[x,], 0.95) })

##plot
pdf(paste(main_path,"/precip_estimates.pdf",sep=""),width=6.0,height=6.0)
ymin<-min(0,0)
ymax<-max(16,16)
yrange<-range(ymin,ymax)
xmin<-min(1,1)  #xmin<-min(1,1)
xmax<-max(100,100)
xrange<-range(xmin,xmax)

date.ticks=c(1,2,10,50,100)
date.labels=c(1,2,10,50,100)

y.ticks=c(0,2,4,6,8,10,12,14,16)
y.labels=c(0,2,4,6,8,10,12,14,16)

####define colour and transparancy
myred <- rgb(1, 102/255, 102/255,0.5)
if("fExtremes" %in% (.packages())){
  detach("package:fExtremes", unload=TRUE)
}
if("evd" %in% (.packages())){
  detach("package:evd", unload=TRUE)
}


plot(1:100,MCrlGTmean, log = "x", xlim = xrange,type="n",xaxt="n",ylim = yrange,yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
axis(1, lwd = 1.5, at=date.ticks, label=date.labels)
axis(2, at = y.ticks, labels = y.labels, lwd = 1.5, las = 1)
polygon(x = c(1:100, 100:1), y = c(upperGT[1:100], rev(lowerGT[1:100])), border = NA , col = "seagreen") #seagreen

par(new=TRUE)
plot(1:100,MCrlMDRmean,log ="x", xlim =xrange,type="n",xaxt="n",ylim = yrange,yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
polygon(x = c(1:100, 100:1), y = c(upperMDR[1:100], rev(lowerMDR[1:100])), border = NA , col = "darkorange") #darkorange

par(new=TRUE)
plot(1:100,MCrlLTmean, log = "x", xlim = xrange,type="n",xaxt="n",ylim = yrange,yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
polygon(x = c(1:100, 100:1), y = c(upperLT[1:100], rev(lowerLT[1:100])), border = NA , col = "red")  


par(new=TRUE)
plot(1:100,MCrlTSmean,log ="x", xlim =xrange,type="n",xaxt="n",ylim = yrange,yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
polygon(x = c(1:100, 100:1), y = c(upperTS[1:100], rev(lowerTS[1:100])), border = NA , col ="blue") #4682B4") #red

par(new=TRUE)
plot(1:100,MCrlNAOmean,log ="x", xlim =xrange,type="n",xaxt="n",ylim = yrange,yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
polygon(x = c(1:100, 100:1), y = c(upperNAO[1:100], rev(lowerNAO[1:100])), border = NA , col = "darkslategray1") ###navy

par(new=TRUE)
plot(1:100,MCrlstatmean,log ="x", xlim =xrange,type="n",xaxt="n",ylim = yrange,yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
polygon(x = c(1:100, 100:1), y = c(upperstat[1:100], rev(lowerstat[1:100])), border = NA , col = "gray70")#wheat

legend(0.65,16.75,legend=c("a) Ellicott City, Maryland"),bty='n',col=1,cex=1.75)#20.25
lines(1:100,MCrlstatmean[1:100],lty=1,lwd=2.0,col="black")
lines(1:100,MCrlGTmean[1:100],lty=2,lwd=2.0,col="seagreen")#seagreen
lines(1:100,MCrlLTmean[1:100],lty=2,lwd=2.0,col="red") #red
lines(1:100,MCrlMDRmean[1:100],lty=2,lwd=2.0,col="darkorange")
lines(1:100,MCrlTSmean[1:100],lty=2,lwd=2.0,col="blue")
lines(1:100,MCrlNAOmean[1:100],lty=2,lwd=2.0,col="darkslategray1")
#####plot the annual maxima observation
par(new=TRUE)
points(1/(1-(1:length(datset))/(length(datset) + 1)), sort(datset), cex =2.5, pch = 20, bg = "black")
par(new=TRUE)
plot(c(2.1,5,10,25,50,98),c(3,4,5,6,8,9),log='x',cex =2.5, pch = 20, bg = "brown4",xpd=TRUE,col="brown4",xlim =xrange,xaxt="n",ylim = yrange,yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
arrows(c(2.1,5,10,25,50,98),noaalow, c(2.1,5,10,25,50,98),noaaup,xpd=TRUE,length=0.05,lwd=3,angle=90,code=3,col="brown4",log='x')

legend("topleft",
       c("Annual maxima observations","NOAA Atlas14 estimate","Stationary","Stationary [5-95% credible interval]","Nonstationary with covariates","Global mean temperature","Atlantic MDR temperature","Local temperature","NAO Index","Atlantic tropical cyclone"),
       col = c("black","brown4","black","black",NA,"black","black","black","black","black"),
       pt.bg = c("black","brown4","black","gray40",NA,"seagreen","darkorange","red","darkslategray1","blue"),
       pch = c(20,20,NA,22,NA,22,22,22,22,22),
       lty = c(NA,NA,1,NA,NA,NA,NA,NA,NA,NA),
       lwd = c(NA,NA,2,NA,NA,NA,NA,NA,NA,NA),
       box.lwd=c(NA,NA,NA,2,NA,2,2,2,2,2),
       pt.cex = c(2.5,2.5,2.5,2.5,NA,2.5,2.5,2.5,2.5,2.5),
       seg.len=c(1,1,1,1,1,1,1,1,1,1), 
        cex=1.1,
         bty="n",
       inset = c(0.01, 0.06))

#############################
y.label='Expected intensity in 2018 [mm/hr]'
x.label='Return period [years]'
xx.label='Probability density'
z.label='Expected 100-yr intensity in 2018 [mm/hr]'
b.label='(b)'
par(new=TRUE)
mtext(y.label, side=2, line=-1.80, outer=TRUE, cex=1.3)
mtext(x.label, side=1, line=-3, adj=0.5,outer=TRUE, cex=1.3)
#mtext(xx.label, side=1, line=-1.1, adj=0.95,outer=TRUE, cex=1.2)
#mtext(z.label, side=4, line=-1.0, outer=TRUE, cex=1.2)
#mtext(z.label, side=2, line=-30.25, outer=TRUE, cex=1.2)
#mtext(b.label,side=1,line=-26.5,adj=0.98,outer=TRUE,cex=1.4)
par(cex=1.0)
 if(postscript==TRUE) {
        dev.off()
}

