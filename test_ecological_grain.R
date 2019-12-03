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
d_cov<-0#(d_sd)^2/2
amu<-(-r/2)
asd<-0.1

##################################
#simulate varying spatial grain, no dispersal
##################################
Nlst<-c(1:30)
niter<-1000
estmat<-as.matrix(data.frame(iter=rep(niter, each=length(Nlst)), tsmp=rep(Nlst, niter),
                             N=NA, n=NA,
                             var=NA, f=NA, r=NA, r_naive=NA, d_sd=NA, d_sd_naive=NA,
                             var_sp=NA, cov_sp=NA))

if(FALSE) {
  n<-1
  for(i in 1:niter) {
    m<-1
    while(m==1 || max(abs(dtot[,-c(1:2)]))>(K*3)) {
      dtot<-symdynN(r = r, amu=amu, asd=asd, f=f, d=d,
                    d_sd=d_sd, d_cov=d_cov, N=max(Nlst),
                    sf=sf, tmax=tmax, stochd = TRUE, stocht = TRUE, as.matrix = TRUE, amax = 0)
      m<-m+1
    }
    dtot<-dtot[dtot[,"time"]>20,]
    cvtmp<-cov(dtot[,-c(1:2),drop=FALSE])
    
    for(j in 1:length(Nlst)) {
      smp<-sample((1:max(Nlst))+2, Nlst[j])
      dtmp<-dtot[,c(1, 2,smp),drop=FALSE]
      estmat[n,"N"]<-Nlst[j]
      estmat[n,"n"]<-nrow(dtmp)
      xsum<-rowSums(dtmp[,-c(1:2),drop=FALSE])
      
      varest<-var(xsum) #mean(xsum^2,na.rm=T)
      estmat[n,"var"]<-varest
      estmat[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1])
      
      #get parameters
      x0<-xsum[-nrow(dtmp)]
      dt<-mean(diff(dtmp[,"time"]))
      ndist<-dtmp[,"disturbed"][-1]
      x1<-xsum[-1]
      x12<-x1^2
      
      d_sd0<-sqrt((Nlst[j]^2-1)*d_cov+Nlst[j]*d_sd^2)
      moddat<-data.frame(x12=x12, x0=x0, dt=dt, ndist=ndist)
      mod2<-try(gnls(sqrt(x12)~sqrt(xt2fun(x0, r=exp(r), d=0, d_sd=exp(d_sd), dt, ndist)), data=moddat,
                    start=c(r=log(r), d_sd=log(d_sd0)), weights = varExp(form=~fitted(.))), silent = TRUE)
      
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
        tmplr<-xsum[disttm_forward+1]/xsum[disttm_forward]
        estmat[n,"r_naive"]<-mean(log(tmplr[tmplr>0])/dt)
      }
      
      #by species
      estmat[n,"var_sp"]<-mean(apply(dtmp[,-c(1:2),drop=FALSE], 2, var))
      if(Nlst[j]>1) {
        cvuse<-cvtmp[smp-2,smp-2]
        estmat[n,"cov_sp"]<-mean(cvuse[row(cvuse)!=col(cvuse)])
      }
      n<-n+1
    }
    if(i/10 == floor(i/10)) {
      cat(paste(round(i/niter,2), ";"))
    }
  }                                        c
  
  write.csv(estmat, "datout/estmat_ec.csv", row.names=F)
 } else {
  estmat<-read.csv("datout/estmat_ec.csv")
}

padj<-c(-1, 0.02, 1.4)

pdf("figures/ecological_grain.pdf", width=6, height=8, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(4,2), mar=c(2,4,1,1), oma=c(2,1.5,0,0))
  Asq<-seq(1, max(Nlst), length=1000)
  
  #r
  ps<-is.finite(estmat[,"r"])
  pltqt(estmat[ps,"N"], estmat[ps,"r"], "", r, domod=FALSE, do_N = FALSE, plog = "", xlab = "", ylim=c(0, 1.1))
  mtext(expression(paste(hat(italic(r)))), 2, line=3.2)
  title("a.", line=padj[1], adj=padj[2], cex.main=padj[3])
  
  ps<-is.finite(estmat[,"r_naive"]) & estmat[,"r_naive"]>0
  pltqt(estmat[ps,"N"], estmat[ps,"r_naive"], "", r, domod=FALSE, do_N = FALSE, plog = "", xlab = "", ylim=c(0, 1.1))
  mtext(expression(paste(italic(r))), 2, line=3.2)
  title("b.", line=padj[1], adj=padj[2], cex.main=padj[3])
  
  #d_sd
  sd_N<-data.frame(Asq, sqrt((d_sd)^2*Asq*(1+(Asq-1)*d_cov/(d_sd)^2)))
  ps<-is.finite(estmat[,"r"])
  pltqt(estmat[ps,"N"], sqrt(estmat[ps,"r"]*(2*f*estmat[ps,"var"])), "", sd_N, domod=FALSE, do_N = FALSE, plog = "", xlab = "", ylim=c(0, 1.8))
  title("c.", line=padj[1], adj=padj[2], cex.main=padj[3])
  mtext(expression(paste(sigma, "|", hat(italic(r)))), 2, line=3.2)
  
  vrest<-mean((estmat[ps,"r"]*(2*f*estmat[ps,"var"]))[estmat[ps,"N"]==1])
  cvest<-mean(estmat[estmat[,"N"]==max(Nlst),"cov_sp"]/estmat[estmat[,"N"]==max(Nlst),"var_sp"],na.rm=T)*vrest
  lines(Asq, sqrt(vrest*Asq+(Asq^2-Asq)*cvest), col=2, lwd=2, lty=3)
  
  ps<-is.finite(estmat[,"d_sd_naive"])
  pltqt(estmat[ps,"N"], estmat[ps,"d_sd_naive"], "", sd_N, domod=FALSE, do_N = FALSE, plog = "", xlab = "", ylim=c(0, 1.8))
  mtext(expression(paste(sigma)), 2, line=3.2)
  title("d.", line=padj[1], adj=padj[2], cex.main=padj[3])
  lines(Asq, sqrt(vrest*Asq+(Asq^2-Asq)*cvest), col=2, lwd=2, lty=3)
  
  #var
  var_N<-data.frame(Asq, var_approx(r,f,d_sd)*Asq*(1+(Asq-1)*d_cov/(d_sd)^2))
  pltqt(estmat[,"N"], estmat[,"var"], "", truev = var_N, do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "", ylim=c(0,5))
  title("e.", line=padj[1], adj=padj[2], cex.main=padj[3])
  mtext(expression(paste("var(", Sigma, italic(x), ")")), 2, line=3.2)
  
  ps<-which(estmat[,"N"]==max(Nlst))
  lines(Asq, Asq*mean(estmat[ps,"var_sp"],na.rm=T)+(Asq^2-Asq)*mean(estmat[ps,"cov_sp"],na.rm=T), col=2, lwd=2, lty=3)
  
  pltqt(estmat[,"N"], estmat[,"var_sp"], "", truev = var_approx(r,f,d_sd), do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "", modoffset = 1e-1, ylim=c(0,0.5))
  mtext(expression(paste("var(", italic(x), ")")), 2, line=3.2)
  title("f.", line=padj[1], adj=padj[2], cex.main=padj[3])
  
  pltqt(estmat[,"N"], estmat[,"cov_sp"]/estmat[,"var_sp"], "", truev = , do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "", modoffset = 1e-1)
  mtext(expression(paste(italic(rho))), 2, line=3.2)
  title("g.", line=padj[1], adj=padj[2], cex.main=padj[3])
  abline(h=0, lty=3)
  
  pltqt(estmat[,"N"], estmat[,"cov_sp"]/estmat[,"var_sp"], "", truev = "", do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "", modoffset = 1e-1, ylim=c(-0.08, 0), qtp = c(qnorm(0.25), qnorm(0.75)))
  abline(h=mean(estmat[estmat[,"N"]==max(Nlst),"cov_sp"]/estmat[estmat[,"N"]==max(Nlst),"var_sp"],na.rm=T), lty=3, col=2, lwd=2)
  mtext(expression(paste(italic(rho))), 2, line=3.2)
  title("h.", line=padj[1], adj=padj[2], cex.main=padj[3])
  mtext("number of species", 1, line = 0.5, outer=T)
  abline(h=0, lty=3)
dev.off()
