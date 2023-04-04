#BIC
statgevbic<-function(data,temps,mu0,sigma_0,xi_0){
  bic<--2*newloglik(data,temps,mu0,sigma_0,xi_0)+(log(length(data))*length(c(mu0,sigma_0,xi_0)))
  return(bic)
}

#AIC
statgevaic<-function(data,temps,mu0,sigma_0,xi_0){
  aic<-(2*length(c(mu0,sigma_0,xi_0)))-2*newloglik(data,temps,mu0,sigma_0,xi_0)
  return(aic)
}

library(bda)

# Log-Likelihood Function
newloglik<-function(dataset,tempset1,mu_0,sigma_0,xi_0)
{
  MU<-mu_0
  SIGMA<-exp(sigma_0)
  XI<-xi_0
  m<-min((1+(XI*(dataset-MU)/SIGMA)))
  if((m<0.00001))return((-Inf))
  else if(any(SIGMA<0.00001))return((-Inf)) #See note 3
  else if(any(XI==0)){
    loglik<-sum(log(
      (1/SIGMA)*
        ((exp(-((dataset-MU)/SIGMA)))**(XI+1))*
        exp(-(exp(-((dataset-MU)/SIGMA))))
    )
    )}
  
  else{
    loglik<-sum(log(
      (1/SIGMA)*
        ((1+(XI*((dataset-MU)/SIGMA)))**(-1-(1/XI)))*
        exp(-(1+XI*((dataset-MU)/SIGMA))**(-1/XI))
    ))      
  } #See note 5
  return(loglik)
} #See note 6


#METROPOLIS HASTINGS ALGORITHM
MetNonstatGEVspsl<-function(dataset,tempset1,n,start,errvect,uvect,sdvect){
  # Start the setup
  cand<-matrix(start,1,3,byrow=TRUE); colnames(cand)<-c("mu0","sigma0","xi0")
  finmat<-matrix(start,1,3,byrow=TRUE); colnames(finmat)<-c("mu0","sigma0","xi0")
  aprobmat<-matrix(NA,1,3,byrow=TRUE); colnames(aprobmat)<-c("mu0","sigma0","xi0")
  errsdmat<-matrix(c(errvect,uvect,sdvect),3,3,byrow=TRUE); colnames(errsdmat)<-c("mu0","sigma0","xi0")
  
  # Metropolis Hastings Algorithm (Random Walk)
  for(i in nrow(finmat)+1:n){
    cand<-rbind(cand,rep(NA,ncol(cand)));finmat<-rbind(finmat,rep(NA,ncol(finmat)));aprobmat<-rbind(aprobmat,rep(NA,ncol(aprobmat)))
    #print(i)
    #mu0
    for (j in 1:ncol(cand)){
      cand[i,j]<-finmat[i-1,j]+rnorm(1,0,errsdmat[1,j]) 
#START      
      if(j==1){
        likely<-exp((newloglik(dataset,tempset1,
                               cand[i,1],
                               finmat[i-1,2],
                               finmat[i-1,3]
        )+dnorm(cand[i,j],mean=uvect[j],sd=sdvect[j],log=TRUE))-
          (newloglik(dataset,tempset1,
                     finmat[i-1,1],
                     finmat[i-1,2],
                     finmat[i-1,3]
          )+dnorm(finmat[i-1,j],mean=uvect[j],sd=sdvect[j],log=TRUE)))
        
        aprobmat[i,j]<-min(1,(likely),na.rm = TRUE)
      }
########################################
if(j==2){
  likely<-exp((newloglik(dataset,tempset1,
                         finmat[i,1],
                         cand[i,2],
                         finmat[i-1,3]
  )+dnorm(cand[i,j],mean=uvect[j],sd=sdvect[j],log=TRUE))-
    (newloglik(dataset,tempset1,
               finmat[i,1],
               finmat[i-1,2],
               finmat[i-1,3]
    )+dnorm(finmat[i-1,j],mean=uvect[j],sd=sdvect[j],log=TRUE)))
  
  aprobmat[i,j]<-min(1,(likely),na.rm = TRUE)
}
########################################
if(j==3){
  likely<-exp((newloglik(dataset,tempset1,
                         finmat[i,1],
                         finmat[i,2],
                         cand[i,3]
  )+dnorm(cand[i,j],mean=uvect[j],sd=sdvect[j],log=TRUE))-
    (newloglik(dataset,tempset1,
               finmat[i,1],
               finmat[i,2],
               finmat[i-1,3]
    )+dnorm(finmat[i-1,j],mean=uvect[j],sd=sdvect[j],log=TRUE)))
  
  aprobmat[i,j]<-min(1,(likely),na.rm = TRUE)
}
#END
      
      
      u<-runif(1)
      if(u<aprobmat[i,j]){finmat[i,j]<-cand[i,j]}
      if(u>=aprobmat[i,j]){finmat[i,j]<-finmat[i-1,j]}
      
    }}
  
  return(list(finmat=finmat, cand=cand, aprobmat=aprobmat, errsdmat=errsdmat, dataset=dataset,tempset1=tempset1))
  
}

##########RESUME #METROPOLIS HASTINGS
ResumeMetNonstatGEVspsl<-function(data,n){
  # Start the setup
  cand<-data$cand;finmat<-data$finmat;aprobmat<-data$aprobmat;errsdmat<-data$errsdmat;dataset<-data$dataset;tempset1<-data$tempset1
  # Metropolis Hastings Algorithm (Random Walk)
  for(i in (nrow(finmat)+1):n){
    cand<-rbind(cand,rep(NA,ncol(cand)));finmat<-rbind(finmat,rep(NA,ncol(finmat)));aprobmat<-rbind(aprobmat,rep(NA,ncol(aprobmat)))
    #print(i)
    for (j in 1:ncol(cand)){
      cand[i,j]<-finmat[i-1,j]+rnorm(1,0,errsdmat[1,j]) 
      #START      
      if(j==1){
        likely<-exp((newloglik(dataset,tempset1,
                               cand[i,1],
                               finmat[i-1,2],
                               finmat[i-1,3]
        )+dnorm(cand[i,j],mean=uvect[j],sd=sdvect[j],log=TRUE))-
          (newloglik(dataset,tempset1,
                     finmat[i-1,1],
                     finmat[i-1,2],
                     finmat[i-1,3]
          )+dnorm(finmat[i-1,j],mean=uvect[j],sd=sdvect[j],log=TRUE)))
        
        aprobmat[i,j]<-min(1,(likely),na.rm = TRUE)
      }
      ########################################
      if(j==2){
        likely<-exp((newloglik(dataset,tempset1,
                               finmat[i,1],
                               cand[i,2],
                               finmat[i-1,3]
        )+dnorm(cand[i,j],mean=uvect[j],sd=sdvect[j],log=TRUE))-
          (newloglik(dataset,tempset1,
                     finmat[i,1],
                     finmat[i-1,2],
                     finmat[i-1,3]
          )+dnorm(finmat[i-1,j],mean=uvect[j],sd=sdvect[j],log=TRUE)))
        
        aprobmat[i,j]<-min(1,(likely),na.rm = TRUE)
      }
      ########################################
      if(j==3){
        likely<-exp((newloglik(dataset,tempset1,
                               finmat[i,1],
                               finmat[i,2],
                               cand[i,3]
        )+dnorm(cand[i,j],mean=uvect[j],sd=sdvect[j],log=TRUE))-
          (newloglik(dataset,tempset1,
                     finmat[i,1],
                     finmat[i,2],
                     finmat[i-1,3]
          )+dnorm(finmat[i-1,j],mean=uvect[j],sd=sdvect[j],log=TRUE)))
        
        aprobmat[i,j]<-min(1,(likely),na.rm = TRUE)
      }
      #END    
      u<-runif(1)
      if(u<aprobmat[i,j]){finmat[i,j]<-cand[i,j]}
      if(u>=aprobmat[i,j]){finmat[i,j]<-finmat[i-1,j]}
      
    }}
  
  return(list(finmat=finmat, cand=cand, aprobmat=aprobmat, errsdmat=errsdmat, dataset=dataset,tempset1=tempset1))
  
}


########################TRACE PLOTS
# Plot #1
MULTnsgevplots<-function(finlist,rm.burn=FALSE,burn=1000){
  finmu0<-finlist$finmat[,1];  finsigma0<-(finlist$finmat[,2]);  finxi0<-finlist$finmat[,3]
  
  if(rm.burn==FALSE){
    par(mfrow=c(1,3), mar=c(2,2,2,2))
    plot(ts(finmu0),ylab="value",xlab="Iteration",main="Markov Chain Results-Mu0")
    abline(v=burn,col="red")
    plot(ts((finsigma0)),ylab="value",xlab="Iteration",main="Markov Chain Results-Sigma0")
    abline(v=burn,col="red")
    plot(ts(finxi0),ylab="value",xlab="Iteration",main="Markov Chain Results-Xi0")
    abline(v=burn,col="red")

    
  }
  else{
    par(mfrow=c(1,3), mar=c(2,2,2,2))
    mu0.burn<-finmu0[burn:length(finmu0)]
    sigma0.burn<-finsigma0[burn:length(finsigma0)]
    xi0.burn<-finxi0[burn:length(finxi0)]
        
    plot(ts(mu0.burn),ylab="value",xlab="Iteration",main="Markov Chain Results-Mu0")
    plot(ts((sigma0.burn)),ylab="value",xlab="Iteration",main="Markov Chain Results-Sigma0")
    plot(ts(xi0.burn),ylab="value",xlab="Iteration",main="Markov Chain Results-Xi0")
  
  }
}

# Rejection Rate
MULTncrej<-function(finlist,burn=1000){
  finmu0<-finlist$finmat[,1];  finsigma0<-(finlist$finmat[,2]);  finxi0<-finlist$finmat[,3]
  
  aprobmu0.burn<-finlist$aprobmat[-(1:burn),1]
  aprobsigma0.burn<-finlist$aprobmat[-(1:burn),2]
  aprobxi0.burn<-finlist$aprobmat[-(1:burn),3]
  
  rejprob<-data.frame(parameter=c( "mu0","sigma0","xi0"),
                      mean_acceptance_rate=c(
                        round(mean(aprobmu0.burn,na.rm=TRUE),4),round(mean(aprobsigma0.burn,na.rm=TRUE),4),round(mean(aprobxi0.burn,na.rm=TRUE),4)
                      ))
  return(rejprob)
}

#############HPD Calculation
hpd <- function(samp,p=0.05)
{
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1:5],]
  return(hpd)
}

###############Return Level
# Return Level Function for GEV
benreturn<-function(t,mu,sigma,xi){
  x<-mu - (sigma/xi)*(1-((-log(1-(1/t)))^(-xi)))
  return(x)
}

benreturnyear<- function(x,mu,sigma,xi){
  t<-(-1)/(exp(-((x-mu)*(xi/sigma)+1)^(-1/xi))-1) 
  return(t)
}

#Density Plots
DensePlot<-function(finlist,burn=1000){
  #######HPD###############
  # Higest Posterior Density. This finds the shortest interval which contains 95% of the MCMC records.
  finmu0<-finlist$finmat[,1];  finsigma0<-(finlist$finmat[,2]);  finxi0<-finlist$finmat[,3]
  
  mu0.burn<-finmu0[burn:length(finmu0)];sigma0.burn<-finsigma0[burn:length(finsigma0)];xi0.burn<-finxi0[burn:length(finxi0)];

  mu0_hpd<-hpd(mu0.burn,p=0.05);mu0_hpd_b<-hpd(mu0.burn,p=0.1)
  sigma0_hpd<-hpd(sigma0.burn,p=0.05);sigma0_hpd_b<-hpd(sigma0.burn,p=0.1)
  xi0_hpd<-hpd(xi0.burn,p=0.05);xi0_hpd_b<-hpd(xi0.burn,p=0.1)
  
  ###############DENSITY PLOTS##########################
  #Calculate the posterio Modes. This are has been sampled at the higest probability. 
  mudat<-data.frame(mu0.x=density(mu0.burn)$x,mu0.y=density(mu0.burn)$y)
  mumode<-mean(mu0.burn)
  sigmadat<-data.frame(sigma0.x=density(sigma0.burn)$x,sigma0.y=density(sigma0.burn)$y)
  sigmamode<-mean(sigma0.burn)
  xidat<-data.frame(xi0.x=density(xi0.burn)$x,xi0.y=density(xi0.burn)$y)
  ximode<-mean(xi0.burn)
  
  ###########TRUE PARAMETERS  
  
  par(mfrow=c(1,3),mar=c(2,2,2,2))
  
  plot(density(mu0.burn),main=expression(paste("Density for sampled ",mu)))
  abline(v = mumode,col="red",lwd=3)
  abline(v = mu0_hpd[1,1],col="gray",lwd=3,lty="dotted")
  abline(v = mu0_hpd[1,2],col="gray",lwd=3,lty="dotted")
  abline(v = mu0_hpd_b[1,1],col="green",lwd=3,lty="dotted")
  abline(v = mu0_hpd_b[1,2],col="green",lwd=3,lty="dotted")
  #legend("topright", 
  #      legend = c("Posterior Mode",
  #                "MLE Estimate",
  #               "90% Credible Int",
  #              "95% Credible Int"), 
  #  fill = c("red","blue","green","gray"), ncol = 1,
  # cex = 1,text.width=5)
  
  
  
  plot(density((sigma0.burn)),main=expression(paste("Density for sampled ",sigma)))
  abline(v = (sigmamode),col="red",lwd=3)
  abline(v = (sigma0_hpd[1,1]),col="gray",lwd=3,lty="dotted")
  abline(v = (sigma0_hpd[1,2]),col="gray",lwd=3,lty="dotted")
  abline(v = (sigma0_hpd_b[1,1]),col="green",lwd=3,lty="dotted")
  abline(v = (sigma0_hpd_b[1,2]),col="green",lwd=3,lty="dotted")
  
  plot(density(xi0.burn),main=expression(paste("Density for sampled ",xi)))
  abline(v = ximode,col="red",lwd=3)
  abline(v = xi0_hpd[1,1],col="gray",lwd=3,lty="dotted")
  abline(v = xi0_hpd[1,2],col="gray",lwd=3,lty="dotted")
  abline(v = xi0_hpd_b[1,1],col="green",lwd=3,lty="dotted")
  abline(v = xi0_hpd_b[1,2],col="green",lwd=3,lty="dotted")

  
}


#############POSTERIOR MODE
posteriormodeGEV<-function(finlist,burn=1000){
  #######HPD###############
  # Higest Posterior Density. This finds the shortest interval which contains 95% of the MCMC records.
  
  finmu0<-finlist$finmat[,1];  finsigma0<-(finlist$finmat[,2]);  finxi0<-finlist$finmat[,3]

  mu0.burn<-finmu0[burn:length(finmu0)];sigma0.burn<-finsigma0[burn:length(finsigma0)];xi0.burn<-finxi0[burn:length(finxi0)];

  ###############DENSITY PLOTS##########################
  #Calculate the posterio Modes. This are has been sampled at the higest probability. 
  mudat<-data.frame(mu0.x=density(mu0.burn)$x,mu0.y=density(mu0.burn)$y)
  mumode<-mean(mu0.burn)
  sigmadat<-data.frame(sigma0.x=density(sigma0.burn)$x,sigma0.y=density(sigma0.burn)$y)
  sigmamode<-mean(sigma0.burn)
  xidat<-data.frame(xi0.x=density(xi0.burn)$x,xi0.y=density(xi0.burn)$y)
  ximode<-mean(xi0.burn)

  
  modeframe=data.frame(pred=c("mu0","sigma0","xi0"),
                       mode=c(mumode,sigmamode,ximode
                       ))
  return(modeframe)
}


#batchmeans
BM_MCMC<-function(finlist,int=1000,burn=0,end=nrow(finlist$finmat),bline=0){
  bmat<-matrix(NA,nrow(finlist$finmat)-burn,ncol(finlist$finmat))
  bmatest<-matrix(NA,nrow(finlist$finmat)-burn,ncol(finlist$finmat))
  source("http://www.stat.psu.edu/~mharan/batchmeans.R")
 #  source("/gpfs/group/kzk10/default/private/susquehanna_hydro/Sanjib/ghub/SourceCode/batchmeans.R")
  # source("personal.psu.edu/muh10/batchmeans.R")
  finmat_bm<-finlist$finmat[burn:end,]
  seqlist<-seq(10,nrow(finmat_bm),by=int)
  for (i in seqlist){
    #print(i)
    for ( j in 1:ncol(finmat_bm)){
      hold<-bm(finmat_bm[(1:i),j])
      bmat[i,j]<-hold$se
      bmatest[i,j]<-hold$est
    }}
 
  bmat<-bmat[!is.na(bmat[,1]),];bmat<-cbind(bmat,seqlist)
  bmatest<-bmatest[!is.na(bmatest[,1]),];bmatest<-cbind(bmatest,seqlist)
  par(mfrow=c(1,3),mar=c(2,2,2,2))
  plot(bmat[,1]~bmat[,4],typ="l",main="BM for mu0");abline(v=bline,col="red")  
  plot(bmat[,2]~bmat[,4],typ="l",main="BM for sigma0");abline(v=bline,col="red")  
  plot(bmat[,3]~bmat[,4],typ="l",main="BM for xi0");abline(v=bline,col="red")  
  return(list(se=bmat,est=bmatest))
  
  }

#batchmeans - Plot Only
plotonlyBM<-function(finlist,burn=1){
  bmat<-finlist$se
  par(mfrow=c(1,3),mar=c(2,2,2,2))
  plot(bmat[,1]~bmat[,4],typ="l",main="BM for mu0")
  abline(v = burn,col="red",lwd=3) 
  plot(bmat[,2]~bmat[,4],typ="l",main="BM for sigma0")
  abline(v = burn,col="red",lwd=3) 
  plot(bmat[,3]~bmat[,4],typ="l",main="BM for xi0")
  abline(v = burn,col="red",lwd=3) 

}

# Sample Means
Mean_MCMC<-function(finlist,int=1000,burn=0,end=nrow(finlist$finmat),bline=0){
  meanmat<-matrix(NA,nrow(finlist$finmat)-burn,ncol(finlist$finmat))
  finmat_mean<-finlist$finmat[burn:end,]
  seqlist<-seq(1,nrow(finmat_mean),by=int)
  for (i in seqlist){
    #print(i)
    for ( j in 1:ncol(finmat_mean)){
      meanmat[i,j]<-mean(finmat_mean[(1:i),j],na.rm=TRUE)
    }}
  meanmat<-meanmat[!is.na(meanmat[,1]),];meanmat<-cbind(meanmat,seqlist)
  par(mfrow=c(1,3),mar=c(2,2,2,2))
  plot(meanmat[,1]~meanmat[,4],typ="l",main="Sample Mean for mu0");abline(v=bline,col="red")  
  plot((meanmat[,2])~meanmat[,4],typ="l",main="Sample Mean for sigma0");abline(v=bline,col="red")  
  plot(meanmat[,3]~meanmat[,4],typ="l",main="Sample Mean for xi0");abline(v=bline,col="red")  
    
  return(meanmat)
  
}

# Sample Means - Plot Only
plotonlyMean<-function(meanmat,burn=1){
  par(mfrow=c(1,3),mar=c(2,2,2,2))
  plot(meanmat[,1]~meanmat[,4],typ="l",main="Sample Mean for mu0")
  abline(v = burn,col="red",lwd=3) 
  plot((meanmat[,2])~meanmat[,4],typ="l",main="Sample Mean for sigma0")
  abline(v = burn,col="red",lwd=3) 
  plot(meanmat[,3]~meanmat[,4],typ="l",main="Sample Mean for xi0")
  abline(v = burn,col="red",lwd=3) 
  }


CredIntervalsGEV<-function(finlist,burn=1000){
  #######HPD###############
  # Higest Posterior Density. This finds the shortest interval which contains 95% of the MCMC records.
  finmu0<-finlist$finmat[,1];  finsigma0<-exp(finlist$finmat[,2]);  finxi0<-finlist$finmat[,3]

  mu0.burn<-finmu0[burn:length(finmu0)];sigma0.burn<-finsigma0[burn:length(finsigma0)];xi0.burn<-finxi0[burn:length(finxi0)];

  mu0_hpd<-hpd(mu0.burn,p=0.05);mu0_hpd_b<-hpd(mu0.burn,p=0.1)
  sigma0_hpd<-hpd(sigma0.burn,p=0.05);sigma0_hpd_b<-hpd(sigma0.burn,p=0.1)
  xi0_hpd<-hpd(xi0.burn,p=0.05);xi0_hpd_b<-hpd(xi0.burn,p=0.1)
  
  ###############DENSITY PLOTS##########################
  #Calculate the posterio Modes. This are has been sampled at the higest probability. 
  mudat<-data.frame(mu0.x=density(mu0.burn)$x,mu0.y=density(mu0.burn)$y)
  mumode<-mean(mu0.burn)
  sigmadat<-data.frame(sigma0.x=density(sigma0.burn)$x,sigma0.y=density(sigma0.burn)$y)
  sigmamode<-mean(sigma0.burn)
  xidat<-data.frame(xi0.x=density(xi0.burn)$x,xi0.y=density(xi0.burn)$y)
  ximode<-mean(xi0.burn)

  
  ####figures
  credints<-data.frame(parameter=c( "mu0","sigma0","xi0"),
                       posterior_mode=c(mumode,sigmamode,ximode),
                       Low_95_CI=c(mu0_hpd[1,1],
                                     sigma0_hpd[1,1],
                                     xi0_hpd[1,1]),
                       
                       High_95_CI=c(mu0_hpd[1,2],
                                      sigma0_hpd[1,2],
                                      xi0_hpd[1,2]))
  return(credints)
  
}


####################BOXPLOT OF ESTIMATES

BoxPlotGEV<-function(finlist,burn=1000){
  #######HPD###############
  # Higest Posterior Density. This finds the shortest interval which contains 95% of the MCMC records.
  dev.off()
  finmu0<-finlist$finmat[,1];  finsigma0<-(finlist$finmat[,2]);  finxi0<-finlist$finmat[,3]
  
  mu0.burn<-finmu0[burn:length(finmu0)];sigma0.burn<-finsigma0[burn:length(finsigma0)];xi0.burn<-finxi0[burn:length(finxi0)];

  boxplot(mu0.burn,sigma0.burn,xi0.burn,names=c("mu0","sigma0","xi0"),
          ylab="Parameter Value",main="Boxplot for Posterior Distributions")
  abline(a = 0,b=0,col="red",lwd=3,lty="dotted") 
  
}

###################SEQUENTIAL
MCMC_breaks_resume<-function(end,intvl,res=FALSE,resdat,datset,tempset1,
                             start_sim,errvect_sim,uvect_sim,sdvect_sim,filetitle){
  if (res==FALSE){
  finrun_sim<-MetNonstatGEVspsl(datset,tempset1,100,start_sim,errvect_sim,uvect_sim,sdvect_sim)
  
  }
  else if(res==TRUE){
    finrun_sim<-resdat
  }
  breaks<-seq(intvl,end,intvl)
  for (i in 1:length(breaks))
  {
    transition_mat<-finrun_sim
    transition_mat$finmat<-finrun_sim$finmat[(nrow(finrun_sim$finmat)-1):nrow(finrun_sim$finmat),]
    transition_mat$cand<-finrun_sim$cand[(nrow(finrun_sim$cand)-1):nrow(finrun_sim$cand),]
    transition_mat$aprobmat<-finrun_sim$aprobmat[(nrow(finrun_sim$aprobmat)-1):nrow(finrun_sim$aprobmat),]
    transition_mat<-ResumeMetNonstatGEVspsl(transition_mat,intvl)
    
    finrun_sim$finmat<-rbind(finrun_sim$finmat,transition_mat$finmat[-(1:2),])
    finrun_sim$cand<-rbind(finrun_sim$cand,transition_mat$cand[-(1:2),])
    finrun_sim$aprobmat<-rbind(finrun_sim$aprobmat,transition_mat$aprobmat[-(1:2),])
    save(finrun_sim,file=paste(filetitle,"_Simulated.Rdata",sep=""))    
  }
  return(finrun_sim)
  }

# thinning
thin_dat<-function(data,thin){
  thin_dat<-data
  thin_dat$finmat<-data$finmat[seq(1,nrow(data$finmat),thin),]
  thin_dat$cand<-data$cand[seq(1,nrow(data$finmat),thin),]
  thin_dat$aprobmat<-data$aprobmat[seq(1,nrow(data$finmat),thin),]
  return(thin_dat)
}

#####################OPTIMIZE VALUES#########################

optimizemat<-function(int,iter,dataset,tempset1,
                      start_sim,errvect_sim,uvect_sim,sdvect_sim){
  
  testparmat<-matrix(ncol=3,nrow=iter,NA)
  testrejmat<-matrix(ncol=3,nrow=iter,NA)
  start_sim<-start_sim+c(rgamma(1,1,1),rgamma(1,1,1),rgamma(1,1,1))
  
  run1<-MetNonstatGEVspsl(dataset,tempset1
                          ,int,start_sim,errvect_sim,uvect_sim,sdvect_sim)
  testrejmat[1,]<-MULTncrej(run1,burn=0)[,2]
  testparmat[1,]<-run1$errsdmat[1,]
  run1<-ResumeMetNonstatGEVspsl(run1,(nrow(run1$finmat)+int))
  testrejmat[2,]<-MULTncrej(run1,burn=((nrow(run1$finmat)-int)))[,2]
  testparmat[2,]<-run1$errsdmat[1,]
  
  i = 3
  while(
    ((any(testrejmat[i-2,]<0.20) | 
      any(testrejmat[i-2,]>0.45) | 
      any(testrejmat[i-1,]<0.20) | 
      any(testrejmat[i-1,]>0.45)) & i<iter)
  ){
    print(i)
    
    run1<-ResumeMetNonstatGEVspsl(run1,(nrow(run1$finmat)+int))
    testrejmat[i,]<-MULTncrej(run1,burn=((nrow(run1$finmat)-int)))[,2]
    testparmat[i,]<-run1$errsdmat[1,]
    
    par(mfrow=c(2,3),mar=c(2,2,2,2))
    apply(run1$finmat,2,plot.ts) # Print trace plots
    #print(apply(run1$finmat,2,effectiveSize)) # Print ESS
    print(testrejmat[i,])         # Print rejection rates
    
    for (j in 1:ncol(testrejmat)){
      if(testrejmat[i,j]>0.45){
        errvect_sim[j]<-errvect_sim[j]*1.5
      }
      if(testrejmat[i,j]<0.20){
        errvect_sim[j]<-errvect_sim[j]*.5
      }
      else{
        errvect_sim[j]<-errvect_sim[j]
      }
    }
    run1$errsdmat[1,]<-errvect_sim
    
    i=i+1
  } #end i while loop
  return(list(rejmat=testrejmat,propmat=testparmat,run1=run1))
}
#########################################################################
### Inverse Gamma Priors
rinvgamma <- function (n, shape, scale = 1) {
  return(1/rgamma(n = n, shape = shape, rate = scale))
}

dinvgamma <- function (x, shape, scale = 1) {
  alpha <- shape
  beta <- scale
  
  
  if(x>0){ 
    log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
    return(log.density)
  }
  else {return(-Inf)}
}



