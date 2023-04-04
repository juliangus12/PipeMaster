#load library
library(Hmisc)
library(gplots)
library(plotrix)
library(ggplot2)
library(Hmisc)
library(gplots)
library(plotrix)
library(ggplot2)
library(evir)
library(plyr)
main_path=getwd()

myreturnlevel <- function(t,mu,sigma,xi){
  x<-qgev(p=1-(1/t),xi=xi,mu=mu,sigma=sigma)
  return(x)}
#--------------------------------------------------------------
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#
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
#####stationary
rtnprd<-load(paste(main_path,"/GTRESULT/stat_widenorm_rtnlevel.RData",sep=''))
retint21<-retint
lower21 <- quantile(retint21,0.05)
median21<-median(retint21)
mean21<-mean(retint21)
upper21 <- quantile(retint21,0.95)
low21<-quantile(retint21,0.25)
upp21<-quantile(retint21,0.75)

####MDR temp
rtnprd<-load(paste(main_path,"/MDRRESULT/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint22<-myreturnlevel(100, paramtime1, (paramtime2), paramtime3)
lower22 <- quantile(retint22,0.05)
median22<-median(retint22)
mean22<-mean(retint22)
upper22 <- quantile(retint22,0.95)
low22<-quantile(retint22,0.25)
upp22<-quantile(retint22,0.75)

###GFDLesm2m_WRF
rtnprd<-load(paste(main_path,"/GFDLesm2m_WRF/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint23<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
lower23 <- quantile(retint23,0.05)
median23<-median(retint23)
mean23<-mean(retint23)
upper23 <- quantile(retint23,0.95)
low23<-quantile(retint23,0.25)
upp23<-quantile(retint23,0.75)


###MPIESMLR_RegCM4
rtnprd<-load(paste(main_path,"/MPIESMLR_RegCM4/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint24<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
lower24 <- quantile(retint24,0.05)
median24<-median(retint24)
mean24<-mean(retint24)
upper24 <- quantile(retint24,0.95)
low24<-quantile(retint24,0.25)
upp24<-quantile(retint24,0.75)

###MPIESMLR_WRF
rtnprd<-load(paste(main_path,"/MPIESMLR_WRF/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint25<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
lower25 <- quantile(retint25,0.05)
median25<-median(retint25)
mean25<-mean(retint25)
upper25 <- quantile(retint25,0.95)
low25<-quantile(retint25,0.25)
upp25<-quantile(retint25,0.75)

####GFDLESM2M
rtnprd<-load(paste(main_path,"/GFDLESM2M/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint26<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
lower26 <- quantile(retint26,0.05)
median26<-median(retint26)
mean26<-mean(retint26)
upper26 <- quantile(retint26,0.95)
low26<-quantile(retint26,0.25)
upp26<-quantile(retint26,0.75)

##CCSM4
rtnprd<-load(paste(main_path,"/CCSM4/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint27<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
lower27 <- quantile(retint27,0.05)
median27<-median(retint27)
mean27<-mean(retint27)
upper27 <- quantile(retint27,0.95)
low27<-quantile(retint27,0.25)
upp27<-quantile(retint27,0.75)

###CanESM2
rtnprd<-load(paste(main_path,"/CanESM2/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint28<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
lower28 <- quantile(retint28,0.05)
median28<-median(retint28)
mean28<-mean(retint28)
upper28 <- quantile(retint28,0.95)
low28<-quantile(retint28,0.25)
upp28<-quantile(retint28,0.75)

###INMCM4
rtnprd<-load(paste(main_path,"/INMCM4/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint29<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
lower29 <- quantile(retint29,0.05)
median29<-median(retint29)
mean29<-mean(retint29)
upper29 <- quantile(retint29,0.95)
low29<-quantile(retint29,0.25)
upp29<-quantile(retint29,0.75)

####MIROCESM
rtnprd<-load(paste(main_path,"/MIROCESM/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint30<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
lower30 <- quantile(retint30,0.05)
median30<-median(retint30)
mean30<-mean(retint30)
upper30 <- quantile(retint30,0.95)
low30<-quantile(retint30,0.25)
upp30<-quantile(retint30,0.75)

###MIROC5
rtnprd<-load(paste(main_path,"/MIROC5/nonstat_param2090.RData",sep=''))
paramtime1<-muproj_chain
paramtime2<-sigma_chain
paramtime3<-xi_chain
retint31<-myreturnlevel(100, paramtime1, paramtime2, paramtime3)
lower31 <- quantile(retint31,0.05)
median31<-median(retint31)
mean31<-mean(retint31)
upper31 <- quantile(retint31,0.95)
low31<-quantile(retint31,0.25)
upp31<-quantile(retint31,0.75)

###
pdf(paste(main_path,"/precip_projections.pdf",sep=""),width=9.0,height=6.0)
####define colour and transparancy
myred <- rgb(1, 102/255, 102/255,0.5)
if("fExtremes" %in% (.packages())){
  detach("package:fExtremes", unload=TRUE)
}
if("evd" %in% (.packages())){
  detach("package:evd", unload=TRUE)
}

xmin<-min(-0.13,-0.13)
xmax<-max(0.65,0.65)
xrange<-range(xmin,xmax)
ymin<-min(4,4)
ymax<-max(20,20)
yrange<-range(ymin,ymax)
date.ticks<-c(4,6,8,10,12,14,16,18,20)
dates<-c(4,6,8,10,12,14,16,18,20)
y.ticks<-c(NA,0,0.45)
ylabels<-c(NA,0,0.45)
##
ret21=density(retint21)
ret22=density(retint22)
ret23=density(retint23)
ret24=density(retint24)
ret25=density(retint25)
ret26=density(retint26)
ret27=density(retint27)
ret28=density(retint28)
ret29=density(retint29)
ret30=density(retint30)
ret31=density(retint31)
plot(ret23$x,ret23$y,col="tomato2",type='l',lty=1,lwd=3.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")
par(new=TRUE)
plot(ret24$x,ret24$y,col="violetred1",type='l',lty=1,lwd=3.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")
par(new=TRUE)
plot(ret25$x,ret25$y,col="sienna",type='l',lty=1,lwd=3.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")
par(new=TRUE)

plot(ret26$x,ret26$y,col="dodgerblue",type='l',lty=1,lwd=3.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")
par(new=TRUE)
plot(ret27$x,ret27$y,col="dodgerblue4",type='l',lty=1,lwd=3.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")
par(new=TRUE)
plot(ret28$x,ret28$y,col="cornflowerblue",type='l',lty=1,lwd=3.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")

par(new=TRUE)
plot(ret29$x,ret29$y,col="lightskyblue",type='l',lty=1,lwd=3.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")

par(new=TRUE)
plot(ret30$x,ret30$y,col="darkslateblue",type='l',lty=1,lwd=3.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")

par(new=TRUE)
plot(ret31$x,ret31$y,col="darkblue",type='l',lty=1,lwd=3.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")

###############
par(new=TRUE)
plot(ret21$x,ret21$y,col="black",type='l',lty=1,lwd=5.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',cex.main=1.0,xaxs="i",yaxs="i")
axis(1, lwd = 1.5, at=date.ticks, label=dates,las=1,cex.axis=1.2)
axis(2, at = y.ticks, labels = ylabels, lwd = 1.5, las = 1,cex.axis=1.2)
par(new=TRUE)
plot(ret22$x,ret22$y,col="black",type='l',lty=2,lwd=5.0,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")
################

par(new=TRUE)
yy1=c(-0.015,-0.010,-0.010,-0.015)
yy2=c(-0.021,-0.025,-0.025,-0.021)
yy3=c(-0.031,-0.035,-0.035,-0.031)
yy4=c(-0.041,-0.045,-0.045,-0.041)
yy5=c(-0.051,-0.055,-0.055,-0.051)
yy6=c(-0.061,-0.065,-0.065,-0.061)
yy7=c(-0.071,-0.075,-0.075,-0.071)
yy8=c(-0.081,-0.085,-0.085,-0.081)
yy9=c(-0.091,-0.095,-0.095,-0.091)
yy10=c(-0.101,-0.105,-0.105,-0.101)
yy11=c(-0.111,-0.115,-0.115,-0.111)

##################
par(new=TRUE)
xx=c(lower21,lower21,upper21,upper21)
zz=c(low21,low21,upp21,upp21)
bb=c(mean21,mean21,mean21,mean21)
aa=c(median21,median21,median21,median21)

polygon(xx,yy1,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy1,border=NA,col="black")
par(new=TRUE)
lines(aa,yy1,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy1,border=NA,lwd=3,col="gray40")

xx=c(lower22,lower22,upper22,upper22)
zz=c(low22,low22,upp22,upp22)
aa=c(median22,median22,median22,median22)
bb=c(mean22,mean22,mean22,mean22)
polygon(xx,yy2,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy2,border=NA,col="gray20")
par(new=TRUE)
lines(aa,yy2,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy2,border=NA,lwd=3,col="gray40")

xx=c(lower23,lower23,upper23,upper23)
zz=c(low23,low23,upp23,upp23)
aa=c(median23,median23,median23,median23)
bb=c(mean23,mean23,mean23,mean23)
polygon(xx,yy10,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy10,border=NA,col="tomato2")
par(new=TRUE)
lines(aa,yy10,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy10,border=NA,lwd=3,col="gray40")


xx=c(lower24,lower24,upper24,upper24)
zz=c(low24,low24,upp24,upp24)
aa=c(median24,median24,median24,median24)
bb=c(mean24,mean24,mean24,mean24)
polygon(xx,yy11,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy11,border=NA,col="violetred1")
par(new=TRUE)
lines(aa,yy11,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy11,border=NA,lwd=3,col="gray40")

xx=c(lower25,lower25,upper25,upper25)
zz=c(low25,low25,upp25,upp25)
aa=c(median25,median25,median25,median25)
bb=c(mean25,mean25,mean25,mean25)
polygon(xx,yy9,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy9,border=NA,col="sienna")
par(new=TRUE)
lines(aa,yy9,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy9,border=NA,lwd=3,col="gray40")

xx=c(lower26,lower26,upper26,upper26)
zz=c(low26,low26,upp26,upp26)
aa=c(median26,median26,median26,median26)
bb=c(mean26,mean26,mean26,mean26)
polygon(xx,yy8,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy8,border=NA,col="dodgerblue")
par(new=TRUE)
lines(aa,yy8,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy8,border=NA,lwd=3,col="gray40")


xx=c(lower27,lower27,upper27,upper27)
zz=c(low27,low27,upp27,upp27)
aa=c(median27,median27,median27,median27)
bb=c(mean27,mean27,mean27,mean27)
polygon(xx,yy7,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy7,border=NA,col="dodgerblue4")
par(new=TRUE)
lines(aa,yy7,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy7,border=NA,lwd=3,col="gray40")

xx=c(lower28,lower28,upper28,upper28)
zz=c(low28,low28,upp28,upp28)
aa=c(median28,median28,median28,median28)
bb=c(mean28,mean28,mean28,mean28)
polygon(xx,yy6,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy6,border=NA,col="cornflowerblue")
par(new=TRUE)
lines(aa,yy6,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy6,border=NA,lwd=3,col="gray40")

xx=c(lower29,lower29,upper29,upper29)
zz=c(low29,low29,upp29,upp29)
aa=c(median29,median29,median29,median29)
bb=c(mean29,mean29,mean29,mean29)
polygon(xx,yy5,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy5,border=NA,col="lightskyblue")
par(new=TRUE)
lines(aa,yy5,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy5,border=NA,lwd=3,col="gray40")

xx=c(lower30,lower30,upper30,upper30)
zz=c(low30,low30,upp30,upp30)
aa=c(median30,median30,median30,median30)
bb=c(mean30,mean30,mean30,mean30)
polygon(xx,yy4,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy4,border=NA,col="darkslateblue")
par(new=TRUE)
lines(aa,yy4,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy4,border=NA,lwd=3,col="gray40")

xx=c(lower31,lower31,upper31,upper31)
zz=c(low31,low31,upp31,upp31)
aa=c(median31,median31,median31,median31)
bb=c(mean31,mean31,mean31,mean31)
polygon(xx,yy3,border=NA,col="gray40")
par(new=TRUE)
polygon(zz,yy3,border=NA,col="darkblue")
par(new=TRUE)
lines(aa,yy3,border=NA,lwd=3,col="black")
par(new=TRUE)
lines(bb,yy3,border=NA,lwd=3,col="gray40")

par(new=TRUE)
plot(9,-0.13,cex=2.5,pch=21,col="brown4",bg="brown4",xpd=TRUE,ylab=NA,xlab=NA,xlim=yrange,ylim=xrange,xaxt='n',yaxt='n',main='',xaxs="i",yaxs="i")
arrows(8,-0.13,10,-0.13,length=0.05,lwd=4,angle=90,code=3,col="brown4",xpd=TRUE)


legend("topleft",
       c("NOAA ATLAS14 estimate","Fitted GEV with observation","Fitted GEV with covariate"),
       col = c("brown4","black","black"),
       lty = c(NA,1,2),
       lwd = c(NA,4,4),
       pch=c(20,NA,NA),
       pt.bg=c("brown4",NA,NA),
       pt.cex=c(3,NA,NA),
       box.lwd=c(NA,NA,NA),
       cex=1.0,
       seg.len=c(1.5,1.5,1),
       inset = c(0.01,0.06))



legend(10.25,0.602,
       c("Global climate models","(Dynamical downscaling)","GFDL-ESM2M.WRF","MPI-ESM-LR.RegCM4","MPI-ESM-LR.WRF"),
       col = c(NA,NA,"tomato2","violetred1","sienna"),
       lty = c(NA,NA,1,1,1),
       lwd = c(NA,NA,3,3,3),
       pch=c(NA,NA,NA,NA,NA),
       pt.bg=c(NA,NA,NA,NA,NA),
       pt.cex=c(NA,NA,NA,NA,NA),
       box.lwd=c(NA,NA,NA,NA,NA),
       cex=1.0,
       seg.len=c(1,1,1,1,1),
       inset = c(0.01, 0.06))


legend("topright",
       c("Global climate models","(Statistical downscaling)","CanESM2","GFDL-ESM2M","CCSM4","INMCM4","MICROCESM","MICRO5"),
       col = c(NA,NA,"cornflowerblue","dodgerblue","dodgerblue4","lightskyblue","darkslateblue","darkblue"),
       lty = c(NA,NA,1,1,1,1,1,1),
       lwd = c(NA,NA,3,3,3,3,3,3),
       pch=c(NA,NA,NA,NA,NA,NA,NA,NA),
       pt.bg=c(NA,NA,NA,NA,NA,NA,NA,NA),
       pt.cex=c(NA,NA,NA,NA,NA,NA,NA,NA),
       box.lwd=c(NA,NA,NA,NA,NA,NA,NA,NA),
       cex=1.0,
       seg.len=c(1,1,1,1,1,1,1),
       inset = c(0.01, 0.06))
legend(3.25,0.68,legend=c("b) Ellicott City, Maryland"),bty='n',col=1,cex=1.75)
#############################
#legend(3.5,-0.05,legend=c("(c)"),bty='n',col=1,cex=1.5)
xx.label='Expected 100-yr rainfall intensity in 2090 [mm/hr]'
#x.label='Return period [years]'
y.label='Probability density'
#z.label='Expected 100-yr intensity in 2018 [mm/hr]'
#b.label='(b)'
par(new=TRUE)
mtext(y.label, side=2, line=-2.0, adj=0.45,outer=TRUE, cex=1.3)
#mtext(x.label, side=1, line=-1.1, adj=0.35,outer=TRUE, cex=1.2)
mtext(xx.label, side=1, line=-2.5, adj=0.45,outer=TRUE, cex=1.3)
#mtext(z.label, side=4, line=-1.0, outer=TRUE, cex=1.2)
#mtext(z.label, side=2, line=-30.25, outer=TRUE, cex=1.2)
#mtext(b.label,side=1,line=-26.5,adj=0.745,outer=TRUE,cex=1.4)
par(cex=1.0)
 if(postscript==TRUE) {
        dev.off()
}

