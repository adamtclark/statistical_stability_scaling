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
sf<-1 #sampling frequency
tmax<-120 #maximum time
d_cov<-d_cov0<-(d_sd)^2/2

##################################
#simulate varying spatial grain, no dispersal
##################################
#set parameters
Alst<-1:30
A<-c(2)
Ilst<-seq(0, 2, by=0.1)
d_cov<-0
niter<-1000

if(FALSE) {
  Alst<-1:30
  niter<-1000
  estmat<-as.matrix(data.frame(iter=rep(niter, each=length(Alst)), tsmp=rep(Alst, niter),
                     N=NA, n=NA,
                     var=NA, f=NA, r=NA, r_naive=NA, d_sd=NA, d_sd_naive=NA))
  
  n<-1
  for(i in 1:niter) {
    dtot<-symdynN(r = r, amu=0, asd=0, f=f, d=d,
                  d_sd=d_sd, d_cov=d_cov, N=max(Alst),
                  sf=sf, tmax=tmax, stochd = TRUE, stocht = TRUE, as.matrix = TRUE, Ifrac=0, dffun = df_col)
    dtot<-dtot[dtot[,"time"]>20,]
    
    for(j in 1:length(Alst)) {
      dtmp<-dtot[,c(1, 2, sample((1:max(Alst))+2, Alst[j])),drop=FALSE]
      estmat[n,"N"]<-Alst[j]
      estmat[n,"n"]<-nrow(dtmp)
      xsum<-rowSums(dtmp[,-c(1:2),drop=FALSE])
      
      varest<-mean(xsum^2,na.rm=T)
      estmat[n,"var"]<-varest
      estmat[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1])
  
      #get parameters
      x0<-xsum[-nrow(dtmp)]
      dt<-mean(diff(dtmp[,"time"]))
      ndist<-dtmp[,"disturbed"][-1]
      x1<-xsum[-1]
      x12<-x1^2
      
      d_sd0<-sqrt((Alst[j]^2-1)*d_cov+Alst[j]*d_sd^2)
      moddat<-data.frame(x12=x12, x0=x0, dt=dt, ndist=ndist)
      mod2<-try(gnls(sqrt(x12)~sqrt(xt2fun(x0, r=exp(r), d=0, d_sd=exp(d_sd), dt, ndist)), data=moddat,
                     start=c(r=log(r+runif(1, min=-0.1, max = 0.1)), d_sd=log(d_sd0)), weights = varExp(form=~fitted(.))), silent = TRUE)
      
      if(!is.character(mod2) & !is.null(mod2)) {
        estmat[n,"d_sd"]<-exp(unname(coef(mod2)["d_sd"]))
        estmat[n,"r"]<-exp(unname(coef(mod2)["r"]))
      }
      
      disttm_backward<-which(c(dtmp[,"disturbed"][-1], 0)==1 & dtmp[,"disturbed"]==0)
      if(sum(disttm_backward)>0) {
        estmat[n,"d_sd_naive"]<-sd(xsum[disttm_backward+1]-xsum[disttm_backward])
      }
      
      disttm_forward<-which(c(dtmp[,"disturbed"][-1], 1)==0 & dtmp[,"disturbed"]==1)
      if(sum(disttm_backward)>0) {
        estmat[n,"r_naive"]<-mean(log(xsum[disttm_forward+1]/xsum[disttm_forward])/dt)
      }
      n<-n+1
    }
    if(i/10 == floor(i/10)) {
      cat(paste(round(i/niter,2), ";"))
    }
  }
  
  ##################################
  #simulate varying spatial grain, with dispersal
  ##################################
  A<-c(2)
  Ilst<-seq(0, 2, by=0.1)
  d_cov<-0
  niter<-1000
  estmat_disp<-as.matrix(data.frame(iter=rep(niter, each=length(Ilst)), tsmp=rep(Ilst, niter),
                               n=NA, I=NA,
                               var=NA, f=NA, r=NA, r_naive=NA, d_sd=NA, d_sd_naive=NA,
                               var_sp=NA, cor_sp=NA))
  n<-1
  for(i in 1:niter) {
    for(j in 1:length(Ilst)) {
      dtmp<-symdynN(r = r, amu=0, asd=0, f=f, d=d,
                    d_sd=d_sd, d_cov=d_cov, N=A,
                    sf=sf, tmax=tmax, stochd = TRUE, stocht = TRUE, as.matrix = TRUE, Ifrac=Ilst[j], dffun = df_col)
      dtmp<-dtmp[dtmp[,"time"]>20,]
      
      estmat_disp[n,"I"]<-Ilst[j]
      estmat_disp[n,"n"]<-nrow(dtmp)
      xsum<-rowSums(dtmp[,-c(1:2),drop=FALSE])
      
      #total
      varest<-mean(xsum^2,na.rm=T)
      estmat_disp[n,"var"]<-varest
      estmat_disp[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1])
      
      #by species
      estmat_disp[n,"var_sp"]<-mean(apply(dtmp[,-c(1:2),drop=FALSE], 2, function(x) mean(x^2)))
      cvtmp<-cov(dtmp[,-c(1:2),drop=FALSE])
      estmat_disp[n,"cor_sp"]<-mean(cvtmp[row(cvtmp)!=col(cvtmp)])
      
      #total
      x0<-xsum[-nrow(dtmp)]
      dt<-mean(diff(dtmp[,"time"]))
      ndist<-dtmp[,"disturbed"][-1]
      x1<-xsum[-1]
      x12<-x1^2
      
      d_sd0<-sqrt((Alst[j]^2-1)*d_cov+Alst[j]*d_sd^2)
      moddat<-data.frame(x12=x12, x0=x0, dt=dt, ndist=ndist)
      mod2<-try(gnls(sqrt(x12)~sqrt(xt2fun(x0, r=exp(r), d=0, d_sd=exp(d_sd), dt, ndist)), data=moddat,
                     start=c(r=log(r+runif(1, min=-0.1, max = 0.1)), d_sd=log(d_sd0)), weights = varExp(form=~fitted(.))), silent = TRUE)
      
      if(!is.character(mod2) & !is.null(mod2)) {
        estmat_disp[n,"d_sd"]<-exp(unname(coef(mod2)["d_sd"]))
        estmat_disp[n,"r"]<-exp(unname(coef(mod2)["r"]))
      }
      
      disttm_backward<-which(c(dtmp[,"disturbed"][-1], 0)==1 & dtmp[,"disturbed"]==0)
      if(sum(disttm_backward)>0) {
        estmat_disp[n,"d_sd_naive"]<-sd(xsum[disttm_backward+1]-xsum[disttm_backward])
      }
      
      disttm_forward<-which(c(dtmp[,"disturbed"][-1], 1)==0 & dtmp[,"disturbed"]==1)
      if(sum(disttm_backward)>0) {
        estmat_disp[n,"r_naive"]<-mean(log(xsum[disttm_forward+1]/xsum[disttm_forward])/dt)
      }
      n<-n+1
    }
    if(i/10 == floor(i/10)) {
      cat(paste(round(i/niter,2), ";"))
    }
  }
  
  write.csv(estmat, "datout/estmat_sp.csv", row.names=F)
  write.csv(estmat_disp, "datout/estmat_sp_disp.csv", row.names=F)
} else {
  estmat<-read.csv("datout/estmat_sp.csv")
  estmat_disp<-read.csv("datout/estmat_sp_disp.csv")
}

pdf("figures/spatial_grain.pdf", width=6, height=6, colormodel = "cmyk", useDingbats = FALSE)
    par(mfcol=c(3,2), mar=c(2,4,1,1), oma=c(2,1.5,0,0))
    #variance and disturbance frequency
    Asq<-seq(1, max(Alst), length=1000)
    var_N<-data.frame(Asq, var_approx(r,f,d_sd)*Asq*(1+(Asq-1)*d_cov0/(d_sd)^2))
    pltqt(estmat[,"N"], estmat[,"var"], "", truev = var_N, do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "spatial grain")
    title("a.", line=-1.2, adj=0.04, cex.main=1.4)
    mtext(expression(paste("var(", Sigma, italic(x), ")")), 2, line=2.7)
    
    #r
    ps<-is.finite(estmat[,"r"])
    pltqt(estmat[ps,"N"], estmat[ps,"r"], "", r, domod=FALSE, do_N = FALSE, plog = "", xlab = "spatial grain", ylim = c(0.9, 1.1), qtp = c(-1, 1))
    title("b.", line=-1.2, adj=0.04, cex.main=1.4)
    mtext(expression(paste(hat(italic(r)))), 2, line=2.7)
    
    #d_sd
    sd_N<-data.frame(Asq, sqrt((d_sd)^2*Asq*(1+(Asq-1)*d_cov0/(d_sd)^2)))
    ps<-is.finite(estmat[,"r"])
    pltqt(estmat[ps,"N"], sqrt(estmat[ps,"r"]*(2*f*estmat[ps,"var"])), "", sd_N, domod=FALSE, do_N = FALSE, plog = "", xlab = "spatial grain")
    title("c.", line=-1.2, adj=0.04, cex.main=1.4)
    mtext(expression(paste(sigma, "|", hat(italic(r)))), 2, line=2.7)
    
    mtext("spatian grain, or number of species", 1, line = 2.5, outer=F)
    
    #variance and disturbance frequency
    var_N<-var_approx(r,f,d_sd)*A*(1+(A-1)*d_cov/(d_sd)^2)
    pltqt(estmat_disp[,"I"], estmat_disp[,"var"], "", truev = var_N, do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "dispersal rate", plotqtl = c(0, 0.975))
    mtext(expression(paste("var(", Sigma, italic(x), ")")), 2, line=2.7)
    title("d.", line=-1.2, adj=0.04, cex.main=1.4)
    
    #r
    ps<-is.finite(estmat_disp[,"r"])
    pltqt(estmat_disp[ps,"I"], estmat_disp[ps,"r"], "", r, domod=FALSE, do_N = FALSE, plog = "", xlab = "dispersal rate", ylim = c(0.9, 1.1), qtp = c(-1, 1))
    mtext(expression(paste(hat(italic(r)))), 2, line=2.7)
    title("e.", line=-1.2, adj=0.04, cex.main=1.4)
    
    #d_sd
    sd_N<-sqrt((d_sd)^2*A*(1+(A-1)*d_cov/(d_sd)^2))
    ps<-is.finite(estmat_disp[,"r"])
    pltqt(estmat_disp[ps,"I"], sqrt(estmat_disp[ps,"r"]*(2*f*estmat_disp[ps,"var"])), "", sd_N, domod=FALSE, do_N = FALSE, plog = "", xlab = "dispersal rate", plotqtl = c(0, 0.975))
    mtext(expression(paste(sigma, "|", hat(italic(r)))), 2, line=2.7)
    title("f.", line=-1.2, adj=0.04, cex.main=1.4)
    mtext("dispersal rate", 1, line = 2.5, outer=F)
dev.off()


pdf("figures/spatial_grain_dispersal_patch_level.pdf", width=3, height=4, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,1), mar=c(2,4,1,1), oma=c(2,1.5,0,0))
  
  pltqt(estmat_disp[,"I"], estmat_disp[,"var_sp"], "", truev = var_approx(r,f,d_sd), do_N = FALSE, domod=FALSE, plog = "", mlog="xy", xlab = "dispersal rate", modoffset = 1e-1)
  mtext(expression(paste("var(", italic(x), ")")), 2, line=3.5)
  title("a.", line=-0.8, adj=0.02, cex.main=1.2)
  
  pltqt(estmat_disp[,"I"], estmat_disp[,"cor_sp"]/estmat_disp[,"var_sp"], "", truev = 0, do_N = FALSE, domod=FALSE, plog = "", mlog="x", xlab = "dispersal rate", modoffset = 1e-1)
  mtext(expression(paste(italic(rho))), 2, line=3.5)
  title("b.", line=-0.8, adj=0.02, cex.main=1.2)
  mtext("dispersal rate", 1, line = 0.5, outer=T, adj=0.65)
dev.off()