library(Hmisc)
library(gplots)
library(plotrix)
library(ggplot2)
library(evir)
library(plyr)

rm(list=ls())

I<-9
C<-0.865
n<-0.013
A1<-0.2
A2<-0.5
A3<-1
S<-0.01
D<-c(18,21,24,27,30,36,42,48)

J1 <- c(18,21,24,27,30,36,42,48)
J2 <- c(18,21,24,27,30,36,42,48)
J3 <- c(18,21,24,27,30,36,42,48)
J4 <- c(18,21,24,27,30,36,42,48)

Pipe <- expand.grid(J1, J2, J3, J4)
colnames(Pipe)[c(1:4)] <- c("J1", "J2", "J3", "J4")
Pipe<- subset(Pipe,Pipe[,2]>=J1)
Pipe<- subset(Pipe,Pipe[,3]>=J2)
Pipe<- subset(Pipe,Pipe[,4]>=J3)

Pipe$cost <- (Pipe$J1/18 + Pipe$J2/18 + Pipe$J3/18 + Pipe$J4/18)

Pipe$Qp1 <-(0.31/n)*((((Pipe$J1)*0.0254)^(8/3))*(S^(1/2)))
Pipe$Qp2 <-(0.31/n)*((((Pipe$J2)*0.0254)^(8/3))*(S^(1/2)))
Pipe$Qp3 <-(0.31/n)*((((Pipe$J3)*0.0254)^(8/3))*(S^(1/2)))
Pipe$Qp4 <-(0.31/n)*((((Pipe$J4)*0.0254)^(8/3))*(S^(1/2)))

Pipe$Qy1<-(0.278*C*I*A1)
Pipe$Qy2<-(0.278*C*I*A2)
Pipe$Qy3<-(0.278*C*I*A3)
Pipe$Qy4<- Pipe$Qy1+Pipe$Qy2+Pipe$Qy3

Pipe$q1<-Pipe$Qp1-Pipe$Qy1
Pipe$q2<-Pipe$Qp2-Pipe$Qy2
Pipe$q3<-Pipe$Qp3-Pipe$Qy3
Pipe$q4<-Pipe$Qp4-Pipe$Qy4


library("dplyr")
stat_1<- Pipe %>% group_by(cost) %>% summarise(n = n()) 
stat_t<- Pipe %>% group_by(cost) %>% add_tally()

stat1<- stat_t %>% group_by(cost) %>% count(TR = q1 < 0)
stat2<- stat_t %>% group_by(cost) %>% count(TR = q2 < 0)
stat3<- stat_t %>% group_by(cost) %>% count(TR = q3 < 0)
stat4<- stat_t %>% group_by(cost) %>% count(TR = q4 < 0)

stat_DF1 <- merge(x=stat1, y=stat_1, by = "cost")
stat_DF2 <- merge(x=stat2, y=stat_1, by = "cost")
stat_DF3 <- merge(x=stat3, y=stat_1, by = "cost")
stat_DF4 <- merge(x=stat4, y=stat_1, by = "cost")

stat_DF1$p<-stat_DF1$p
stat_DF2$p<-stat_DF2$p
stat_DF3$p<-stat_DF3$p
stat_DF4$p<-stat_DF4$p

stat_DF1$p<- stat_DF1$n.x/stat_DF1$n.y
stat_DF2$p<- stat_DF2$n.x/stat_DF2$n.y
stat_DF3$p<- stat_DF3$n.x/stat_DF3$n.y
stat_DF4$p<- stat_DF4$n.x/stat_DF4$n.y

stat_DF1$p[stat_DF1$TR==TRUE]<- 0
stat_DF2$p[stat_DF2$TR==TRUE]<- 0
stat_DF3$p[stat_DF3$TR==TRUE]<- 0
stat_DF4$p[stat_DF4$TR==TRUE]<- 0



set.seed(200)                                            
##
# plot(stat_DF1$cost,stat_DF1$p, main="Cost vs Reliability (RI-10mm/hr)", ylab="Reliability of Pipe System (%)", xlab="Cost Factor")
# points(stat_DF2$cost,stat_DF2$p,col='red')
# points(stat_DF3$cost,stat_DF3$p,col='violet')
# points(stat_DF4$cost,stat_DF4$p,col='green')

##
par(mfrow=c(2,2))
plot(stat_DF1$cost,1-stat_DF1$p, main="Catchment Area Junction 1 (RI-9mm/hr)", ylab="Probability of overflow at J1", xlab="Cost (Function of Combined Pipe Diameters)")
plot(stat_DF2$cost,1-stat_DF2$p, col='red',main="Catchment Area Junction 2 (RI-9mm/hr)", ylab="Probability of overflow at J2", xlab="Cost (Function of Combined Pipe Diameters)")
plot(stat_DF3$cost,1-stat_DF3$p, col='violet',main="Catchment Area Junction 3 (RI-9mm/hr)", ylab="Probability of overflow at J3", xlab="Cost (Function of Combined Pipe Diameters)")
plot(stat_DF4$cost,1-stat_DF4$p, col='green',main="Combined Flow Junction 4  (RI-9mm/hr)", ylab="Probability of overflow at J4", xlab="Cost (Function of Combined Pipe Diameters)")
#-----------------------------------------------------#






