#plots Figure 3 in the main text

#load functions and packages
source("functions.R")
require(mvtnorm)
require(nlme)
require(deSolve)

#parameters
K<-1 #carrying capacity
r<-1 #rate of recovery
d<-(0) #mean size of disturbance (mu in text)
d_sd<-sqrt(0.1) #SD of disturbances (sigma in text)
f<-1 #average time between disturbances (1/lambda in text)
sf<-0.01 #sampling frequency
tmax<-120 #maximum time for simulation

set.seed(1202) #set random seed

##################################
#simulate varying sampling frequency
##################################
tlst<-c(0.01, seq(0.05, 10, by=0.5), 10) #sampling frequencies to test
niter<-1000 #number of iterations
#dataframe for saving results
estmat<-as.matrix(data.frame(iter=rep(niter, each=length(tlst)), tsmp=rep(tlst, niter),
                   n=NA,
                   var=NA, f=NA, r=NA, d_sd=NA, r_naive=NA, d_sd_naive=NA))

#only run if no saved simulation is available
if(sum(grep("estmat_tm.csv", dir("datout/")))==0) {
  n<-1
  for(i in 1:niter) {
    for(j in 1:length(tlst)) {
      dtmp<-symdyn(r, f, d, d_sd, sf=tlst[j], tmax=tmax, stochd = TRUE, stocht = TRUE, as.matrix = TRUE)
      dtmp<-dtmp[dtmp[,"time"]>20,] #throw out first 20 time steps as burn-in
      
      estmat[n,"n"]<-nrow(dtmp) #number of observations
      varest<-mean(dtmp[,"state"]^2) #variance (recall E[x]=0)
      estmat[n,"var"]<-varest
      estmat[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1]) #mean realized time between disturbances
  
      #estimate r as a log ratio of change per time
      #this is the "raw" estimate from the paper
      tmp<-try(exp(mean(log(abs(diff(log(abs(dtmp[,"state"]))))),na.rm=T))/mean(diff(dtmp[,"time"]),na.rm=T))
      if(!is.character(tmp)) {
        estmat[n,"r_naive"]<-tmp
      }
      
      #get corrected estimates
      x0<-dtmp[,"state"][-nrow(dtmp)] #observations at time t
      dt<-mean(diff(dtmp[,"time"])) #time step size in sampling
      ndist<-dtmp[,"disturbed"][-1] #number of disturbances
      x1<-dtmp[,"state"][-1] #observations at time t+1
      x12<-x1^2 #x1 squared
      
      #try fitting r, and d_sdbased on xt2fun - this is Eq. S39 from the supplement
      mod2<-try(nls(sqrt(x12)~sqrt(xt2fun(x0, r=exp(r), d=0, d_sd=exp(d_sd), dt, ndist)),
              start=c(r=log(r), d_sd=log(d_sd))), silent = TRUE)
      
      #if model converges, save paramters
      if(!is.character(mod2) & !is.null(mod2)) {
        estmat[n,"d_sd"]<-exp(unname(coef(mod2)["d_sd"]))
        estmat[n,"r"]<-exp(unname(coef(mod2)["r"]))
      }
      
      #estiamte d_sd by comparing time steps with and without disturbances
      #this is the "raw" estimate from the paper
      disttm_backward<-which(c(dtmp[,"disturbed"][-1], 0)==1 & dtmp[,"disturbed"]==0)
      if(sum(disttm_backward)>0) {
        estmat[n,"d_sd_naive"]<-sd(dtmp[disttm_backward+1,"state"]-dtmp[disttm_backward,"state"])
      }
      n<-n+1
    }
    if(i/10 == floor(i/10)) {
      cat(paste(round(i/niter,2), ";"))
    }
  }
  write.csv(estmat, "datout/estmat_tm.csv", row.names=F)
} else {
  #load saved simulation
  estmat<-read.csv("datout/estmat_tm.csv")
}

#upper quantile for plot margins
pdf("figures/temporal_grain.pdf", width=3, height=6, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(3,1), mar=c(2,4,1,1), oma=c(2,1.5,0,0))
  
  #corrected r
  ps<-is.finite(estmat[,"r"])
  pltqt(estmat[ps,"tsmp"], estmat[ps,"r"], "", r, domod=FALSE, do_N = FALSE, plog = "", ylim = c(0, 3.2))
  mtext(expression(paste("resilience, ", italic(r))), 2, line=3.2)
  title("a.", line=-0.85, adj=0.02, cex.main=1.2)
  abline(h=0, lty=3)
  
  #raw r estimate
  ps<-is.finite(estmat[,"r_naive"])
  addqt(estmat[ps,"tsmp"], estmat[ps,"r_naive"])
  
  #corrected d_sd
  ps<-is.finite(estmat[,"r"])
  #back-calculate from corrected r value
  pltqt(estmat[ps,"tsmp"], sqrt(estmat[ps,"r"]*(2*f*estmat[ps,"var"])), "", d_sd, domod=FALSE, do_N = FALSE, plog = "", mlog="", ylim = c(0, 0.6))
  mtext(expression(paste("resistance"^{-1}, ", ", italic(sigma))), 2, line=3.2)
  title("b.", line=-0.85, adj=0.02, cex.main=1.2)
  abline(h=0, lty=3)
  
  #raw d_sd
  ps<-is.finite(estmat[,"d_sd_naive"])
  addqt(estmat[ps,"tsmp"], estmat[ps,"d_sd_naive"])
  
  #variance
  pltqt(estmat[,"tsmp"], estmat[,"var"], "", var_approx(r,f,d_sd), do_N = FALSE, domod=FALSE, plog = "", mlog="", ylim = c(0, 0.1))
  mtext(expression(paste("invariability"^-1, ", ", "var(", italic(x), ")")), 2, line=3.2)
  title("c.", line=-0.85, adj=0.02, cex.main=1.2)
  abline(h=0, lty=3)
  
  mtext("sampling frequency", 1, line = 0.5, outer=T, adj = 0.65)
dev.off()
