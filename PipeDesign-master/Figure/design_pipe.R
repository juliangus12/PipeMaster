
#####load library
library(RColorBrewer)
library(Hmisc)
library(gplots)
library(plotrix)
library(ggplot2)
main_path=getwd()
infile<-paste(main_path,"/dataset/design_pipe.txt",sep='')
big<-scan(infile,skip=4,list(dia1=0,prob1=0,dia2=0,prob2=0,dia3=0,prob3=0,dia4=0,prob4=0,dia5=0,prob5=0,dia6=0,prob6=0,dia7=0,prob7=0,dia8=0,prob8=0))
dia1<-big$dia1
prob1<-big$prob1
dia2<-big$dia2
prob2<-big$prob2
dia3<-big$dia3
prob3<-big$prob3
dia4<-big$dia4
prob4<-big$prob4
dia5<-big$dia5
prob5<-big$prob5
dia6<-big$dia6
prob6<-big$prob6
dia7<-big$dia7
prob7<-big$prob7
dia8<-big$dia8
prob8<-big$prob8
##save output
outfile<-paste("design_pipe.pdf",sep="")
postscript(outfile,horizontal=FALSE)
par(mfrow=c(1,1))
par(oma=c(26,2.1,2.7,1.75))
par(pin=c(6.2,5.5))
par(cex=0.8)
xmin<-min(0.5,0.5)
xmax<-max(7.5,7.5)
xrange<-range(xmin,xmax)
ymin<-min(0.00008,0.00008)
ymax<-max(1.4,1.4)
yrange<-range(ymin,ymax)
x.ticks<-c(1,2,3,4,5,6,7)
x.labels<-c(21,24,27,30,36,42,48)  
xx.labels<-c(1,1.2,1.4,2.4,3.2,3.7,4.4)
y.ticks<-c(0.0001,0.01,0.10,0.39,1.0)
y.labels<-c('1',0.99,0.90,0.61,0)
yy.ticks<-c(0.76,0.076)
yy.labels<-c(1,0.2)
######
dia<-c(1,2,3,4,5,6,7)
proba<-c(prob2[1],prob3[1],prob4[1],prob5[1],prob6[1],prob7[1],prob8[1])
probb<-c(prob2[2],prob3[2],prob4[2],prob5[2],prob6[2],prob7[2],prob8[2])
probc<-c(prob2[3],prob3[3],prob4[3],prob5[3],prob6[3],prob7[3],prob8[3])
probd<-c(prob2[4],prob3[4],prob4[4],prob5[4],prob6[4],prob7[4],prob8[4])
probe<-c(prob2[5],prob3[5],prob4[5],prob5[5],prob6[5],prob7[5],prob8[5])
probf<-c(prob2[6],prob3[6],prob4[6],prob5[6],prob6[6],prob7[6],prob8[6])
probg<-c(prob2[7],prob3[7],prob4[7],prob5[7],prob6[7],prob7[7],prob8[7])
probh<-c(prob2[8],prob3[8],prob4[8],prob5[8],prob6[8],prob7[8],prob8[8])
probi<-c(prob2[9],prob3[9],prob4[9],prob5[9],prob6[9],prob7[9],prob8[9])
probj<-c(prob2[10],prob3[10],prob4[10],prob5[10],prob6[10],prob7[10],prob8[10])
probk<-c(prob2[1],prob3[11],prob4[11],prob5[11],prob6[11],prob7[11],prob8[11])
####
diariab=c(1-0.5,2,3,4,5,6+0.5)
riab<-c(0.39,0.39,0.39,0.39,0.39,0.39)
plot(diariab,riab,type="l",log='y',lty=1,lwd=2.0,col="darkgreen",xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
xx1=c(0,8,8,0)
yy4=c(0.00008,0.00008,0.076,0.076)
yy3=c(0.076,0.076 ,0.76,0.76)
xx=c(0,8,8,0)
yy=c(0.76,0.76,1.4,1.4)
##polygon
polygon(xx, yy4, border = NA , col = "darkseagreen")
par(new=TRUE)
polygon(xx, yy3, border = NA , col = "khaki")
par(new=TRUE)
polygon(xx, yy, border = NA , col = "lightpink")
par(new=TRUE)
###plot reliability
plot(dia,proba,type='p',log='y',pch=21,cex=3.0,bg="black",col="black",lwd=0.9,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
axis(1,at=x.ticks,labels=x.labels,cex.axis=1.75,las=1)
axis(2,at=y.ticks,labels=y.labels,cex.axis=1.75,las=0)
axis(4,at=yy.ticks,labels=yy.labels,cex.axis=1.75,las=0)
axis(3,at=x.ticks,labels=xx.labels,cex.axis=1.75,las=1)
par(new=TRUE)
plot(dia,probb,type='p',log='y',pch="-",col="black",cex=5.0,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
par(new=TRUE)
plot(dia,probc,type='p',log='y',pch="-",cex=5.0,col="tomato2",lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
par(new=TRUE)
plot(dia,probd,type='p',log='y',pch="-",col="violetred1",cex=5.0,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
par(new=TRUE)
plot(dia,probe,type='p',log='y',pch="-",col="sienna",cex=5.0,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
par(new=TRUE)
plot(dia,probf,type='p',log='y',pch='-',col="gray40",cex=5.0,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
par(new=TRUE)
plot(dia,probg,type='p',log='y',pch='-',col="gray45",cex=5.0,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
par(new=TRUE)
plot(dia,probh,type='p',log='y',pch='-',col="gray50",cex=5.0,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
par(new=TRUE)
plot(dia,probi,type='p',log='y',pch='-',col="gray55",cex=5.0,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
par(new=TRUE)
plot(dia,probj,type='p',log='y',pch='-',col="gray60",cex=5.0,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
par(new=TRUE)
plot(dia,probk,type='p',log='y',pch='-',col="gray65",cex=5.0,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange),xaxs="i",yaxs="i")
points(1,0.42,type='p',pch=25,bg='brown4',col="brown4",cex=3,lty='solid',xpd=TRUE,xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange))
par(new=TRUE)
plot(0.77,1.23,type='p',pch="|",col="brown4",cex=2.0,lwd=4,lty='solid',xaxt="n",yaxt="n",ylab='',xlab='',xlim=xrange,ylim=rev(yrange))
arrows(0.8,1.44,4.0,1.44,code=2,xpd=NA,lwd=5,col="black")
text(4.45,1.425,log='y',"SF=1.4",col="black",ncol=2,cex=1.6)
##add text
par(new=TRUE)
text(1.2,1.17,log='y',"Design diameter",col="brown4",ncol=2,cex=1.6)
text(7,0.972,log='y',"<1:500 AEP",col="darkgreen",ncol=2,cex=1.6)
text(7,1.39,log='y',">1:100 AEP",col="red",ncol=2,cex=1.6)
par(new=TRUE)
#legend
legend("topleft",
       c("Fitted GEV with observation"," ","Observation","","GFDL-ESM2M.WRF","MPI-ESM-LR.RegCM4","MPI-ESM-LR.WRF","","CCSM4","CanESM2","GFDL-ESM2M","INMCM4","MIROC5","MIROCESM"),
     col = c("black",NA,"black",NA,"tomato2","violetred1","sienna",NA,"gray40","gray45","gray50","gray55","gray60","gray65"),
       pt.bg = c("black",NA,"black",NA,"tomato2","violetred1","sienna",NA,"gray40","gray45","gray50","gray55","gray60","gray65"),
       pch = c(21,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
       lty = c(NA,NA,2,NA,2,2,2,NA,2,2,2,2,2,2),
       lwd = c(NA,NA,5,NA,5,5,5,NA,5,5,5,5,5,5),
       box.lwd=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
       cex=1.4,
       x.intersp=-0.2,
       pt.cex = c(2,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
       inset = c(0.01, -0.01))

text(1.75,0.05,log='y',"Fitted GEV with covariate",col="black",ncol=2,cex=1.6)
text(1.70,0.18,log='y',"Regional climate models",col="black",ncol=2,cex=1.6)
text(1.6,0.42,log='y',"Global climate models",col="black",ncol=2,cex=1.6)
par(new=TRUE)
##label
y.label='Hydraulic reliability'
x.label='Pipe diameter [inches]'
z.label='Annual exceedance probability [%]'
xx.label='Cost factor'
#xxx.label='a) Ellicott City, Maryland'
par(new=TRUE)
mtext(y.label, side=2, line=0.4, outer=TRUE, cex=1.75)
mtext(x.label, side=1, line=-0.15, outer=TRUE, cex=1.75)
mtext(z.label, side=4, line=0.3, outer=TRUE, cex=1.75)
mtext(xx.label, side=3, line=-0.8, outer=TRUE, cex=1.75)
#mtext(xxx.label,font=1,side=3,line=1.0,adj=-0.05,outer=TRUE,cex=2.0)
par(cex=1.0)
 if(postscript==TRUE) {
        dev.off()
}


