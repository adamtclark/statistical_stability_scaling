error
rm(list=ls())
setwd("~/Dropbox/Projects/052_GFO_2019/src/simulation_study/")
source("functions.R")
source("~/Dropbox/Rfunctions/logit_funs.R")
require(mvtnorm)
require(nlme)
require(deSolve)

set.seed(1202)

K<-1 #carrying capacity
r<-1 #rate of recovery
d<-(0) #mean size of disturbance
d_sd<-sqrt(0.1) #SD of disturbances
f<-1 #mean frequency of disturbance
sf<-0.01 #sampling frequency
tmax<-120 #maximum time


##################################
#simulate varying temporal grain
##################################
tlst<-c(0.01, seq(0.05, 10, by=0.5), 10)
niter<-1000
estmat<-as.matrix(data.frame(iter=rep(niter, each=length(tlst)), tsmp=rep(tlst, niter),
                   n=NA, #neff=NA,
                   var=NA, f=NA, r=NA, d_sd=NA, r_naive=NA, d_sd_naive=NA))

if(FALSE) {
  n<-1
  for(i in 1:niter) {
    for(j in 1:length(tlst)) {
      dtmp<-symdyn(r, f, d, d_sd, sf=tlst[j], tmax=tmax, stochd = TRUE, stocht = TRUE, as.matrix = TRUE)
      dtmp<-dtmp[dtmp[,"time"]>20,]
      
      estmat[n,"n"]<-nrow(dtmp)
      varest<-mean(dtmp[,"state"]^2)
      estmat[n,"var"]<-varest
      estmat[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1])
  
      tmp<-try(exp(mean(log(abs(diff(log(abs(dtmp[,"state"]))))),na.rm=T))/mean(diff(dtmp[,"time"]),na.rm=T))
      if(!is.character(tmp)) {
        estmat[n,"r_naive"]<-tmp
      }
      
      #get parameters
      x0<-dtmp[,"state"][-nrow(dtmp)]
      dt<-mean(diff(dtmp[,"time"]))
      ndist<-dtmp[,"disturbed"][-1]
      x1<-dtmp[,"state"][-1]
      x12<-x1^2
      
      mod2<-try(nls(sqrt(x12)~sqrt(xt2fun(x0, r=exp(r), d=0, d_sd=exp(d_sd), dt, ndist)),
              start=c(r=log(r), d_sd=log(d_sd))), silent = TRUE)
      
      if(!is.character(mod2) & !is.null(mod2)) {
        estmat[n,"d_sd"]<-exp(unname(coef(mod2)["d_sd"]))
        estmat[n,"r"]<-exp(unname(coef(mod2)["r"]))
      }
      
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
  estmat<-read.csv("datout/estmat_tm.csv")
}

qtlu<-0.975

pdf("figures/temporal_grain.pdf", width=3, height=6, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(3,1), mar=c(2,4,1,1), oma=c(2,1.5,0,0))
  #r
  ps<-is.finite(estmat[,"r"])
  pltqt(estmat[ps,"tsmp"], estmat[ps,"r"], "", r, domod=FALSE, do_N = FALSE, plog = "", ylim = c(0, quantile(estmat[ps,"r"], qtlu)))
  mtext(expression(paste("resilience, ", italic(r))), 2, line=3.2)
  title("a.", line=-0.85, adj=0.02, cex.main=1.2)
  abline(h=0, lty=3)
  
  #r naive
  ps<-is.finite(estmat[,"r_naive"])
  addqt(estmat[ps,"tsmp"], estmat[ps,"r_naive"])
  
  # sd
  ps<-is.finite(estmat[,"r"])
  pltqt(estmat[ps,"tsmp"], sqrt(estmat[ps,"r"]*(2*f*estmat[ps,"var"])), "", d_sd, domod=FALSE, do_N = FALSE, plog = "", mlog="", ylim = c(0, quantile(sqrt(estmat[ps,"r"]*(2*f*estmat[ps,"var"])), qtlu)))
  mtext(expression(paste("resistance"^{-1}, ", ", italic(sigma))), 2, line=3.2)
  title("b.", line=-0.85, adj=0.02, cex.main=1.2)
  abline(h=0, lty=3)
  
  #naive sd
  ps<-is.finite(estmat[,"d_sd_naive"])
  addqt(estmat[ps,"tsmp"], estmat[ps,"d_sd_naive"])
  
  #variance and disturbance frequency
  pltqt(estmat[,"tsmp"], estmat[,"var"], "", var_approx(r,f,d_sd), do_N = FALSE, domod=FALSE, plog = "", mlog="", ylim = c(0, quantile(estmat[,"var"], qtlu)))
  mtext(expression(paste("variability, ", "var(", italic(x), ")")), 2, line=3.2)
  title("c.", line=-0.85, adj=0.02, cex.main=1.2)
  abline(h=0, lty=3)
  
  mtext("temporal frequency", 1, line = 0.5, outer=T, adj = 0.65)
dev.off()
