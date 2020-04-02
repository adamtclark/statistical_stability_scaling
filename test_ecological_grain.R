#plots Figure 5 in the main text
#rm(list=ls())

#load functions and packages
source("functions.R")
source("~/Dropbox/Rfunctions/logit_funs.R")
require(mvtnorm)
require(nlme)
require(deSolve)
collst_attributes<-c("firebrick3", "darkgoldenrod3", "dodgerblue2", "mediumpurple2", "olivedrab4", "springgreen4")
#resistance; resilience; invariability; frequency; persistence; robustness

set.seed(1202) #set random seed

K<-1 #carrying capacity
r<-1 #rate of recovery
d<-(0) #mean size of disturbance (mu in text)
d_sd<-sqrt(0.1) #SD of disturbances (sigma in text)
f<-1 #average time between disturbances (1/lambda in text)
sf<-0.1 #sampling interval
tmax<-120 #maximum time for simulation
d_cov<-0 #covariance in disturbances among species
amu<-(-r/2) #average interaction coefficient
asd<-0.1 #standard deviation of interaction coefficient

##################################
#simulate varying ecological grain, no dispersal
##################################
Nlst<-c(1:30) #species numbers to test
niter<-1000 #number of iterations
estmat<-as.matrix(data.frame(iter=rep(niter, each=length(Nlst)), tsmp=rep(Nlst, niter),
                             N=NA, n=NA,
                             var=NA, f=NA, r=NA, r_naive=NA, r_naive_median=NA, r_naive_median_bounded=NA,
                             d_sd=NA, d_sd_naive=NA,
                             var_sp=NA, cov_sp=NA,
                             var_dist=NA, cov_dist=NA,
                             d_sd_true=NA,
                             sumcovmat_sp=NA, sumcovmat_dist=NA))

#only run if no saved simulation is available
if(sum(grep("estmat_ec.csv", dir("datout/")))==0) {
  #simulate samplings of different numbers of species
  n<-1
  for(i in 1:niter) {
    m<-1
    #exclude runs where system "blows up"
    while(m==1 || max(abs(dtot[,-c(1:2)]))>(K*3)) {
      tmpout<-symdynN(r = r, amu=amu, asd=asd, f=f, d=d,
                    d_sd=d_sd, d_cov=d_cov, N=max(Nlst),
                    sf=sf, tmax=tmax, stochd = TRUE, stocht = TRUE, fullout = TRUE, amax = 0)
      dtot<-as.matrix(tmpout$datout)
      m<-m+1
    }
    dtot<-dtot[dtot[,"time"]>20,] #throw out burn-in time
    cvtmp<-cov(dtot[,-c(1:2),drop=FALSE]) #observed covariance of species abundances
    cvdist<-cov(tmpout$dquant) #observed covariance of disturbances
    
    for(j in 1:length(Nlst)) {
      #sample subset of species
      smp<-sample((1:max(Nlst))+2, Nlst[j])
      dtmp<-dtot[,c(1, 2,smp),drop=FALSE]
      estmat[n,"N"]<-Nlst[j]
      estmat[n,"n"]<-nrow(dtmp)
      #observation is summed abundance
      xsum<-rowSums(dtmp[,-c(1:2),drop=FALSE])
      
      #observed variance
      varest<-var(xsum)
      estmat[n,"var"]<-varest
      #sampling frequency
      estmat[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1])
      estmat[n,"d_sd_true"]<-sd(rowSums(tmpout$dquant[,smp-2,drop=FALSE])) #achieved disturbance standard deviations
      
      #estimate parameters
      x0<-xsum[-nrow(dtmp)]
      dt<-mean(diff(dtmp[,"time"]))
      ndist<-dtmp[,"disturbed"][-1]
      x1<-xsum[-1]
      x12<-x1^2
      
      #raw estimate of r
      disttm_forward<-which(c(dtmp[,"disturbed"][-1], 1)==0 & dtmp[,"disturbed"]==1)
      if(sum(disttm_forward)>0) {
        ps<-which(sign(xsum[disttm_forward+1]) == sign(xsum[disttm_forward+1])) &
          (abs(xsum[disttm_forward])>0.1)
        tmplr<-abs(xsum[disttm_forward+1][ps])/abs(xsum[disttm_forward][ps])
        estmat[n,"r_naive"]<-(-mean(log(tmplr)/dt))
        estmat[n,"r_naive_median_bounded"]<-(-median(log(tmplr)/dt))
      }
      #median
      tmplr_tot<-(xsum[disttm_forward+1])/(xsum[disttm_forward])
      tmplr_tot[tmplr_tot<0]<-NA #remove sign changes
      rest<-(-(log(tmplr_tot)/dt))
      estmat[n,"r_naive_median"]<-median(rest, na.rm=T)
      
      #species-level variance and between-species covariance
      estmat[n,"var_sp"]<-mean(diag(cvtmp[smp-2,smp-2]))
      if(Nlst[j]>1) {
        cvuse<-cvtmp[smp-2,smp-2]
        estmat[n,"cov_sp"]<-mean(cvuse[row(cvuse)!=col(cvuse)])
      }
      
      #variance and between-species covariance of disturbances
      estmat[n,"var_dist"]<-mean(diag(cvdist[smp-2,smp-2]))
      if(Nlst[j]>1) {
        cvuse<-cvdist[smp-2,smp-2]
        estmat[n,"cov_dist"]<-mean(cvuse[row(cvuse)!=col(cvuse)])
      }
      
      #summed covariance matrices for calculating r
      estmat[n,"sumcovmat_sp"]<-sum(cvtmp[smp-2,smp-2])
      estmat[n,"sumcovmat_dist"]<-sum(cvdist[smp-2,smp-2])
      
      n<-n+1
    }
    if(i/10 == floor(i/10)) {
      cat(paste(round(i/niter,2), ";"))
    }
  }                                        
  
  write.csv(estmat, "datout/estmat_ec.csv", row.names=F)
} else {
  estmat<-read.csv("datout/estmat_ec.csv")
}

#get r vs. time:
set.seed(1202)
ntm<-10
dtp<-0.1

if(sum(grep("estmat_ec_r.csv", dir("datout/")))==0) {
  #calculate r for different sampling frequencies
  rmat<-matrix(nrow=niter, ncol=ntm/dtp)
  for(ii in 1:niter) {
    m<-1
    while(m==1 || max(abs(dtot[,-c(1:2)]))>(K*3)) {
      #burn-in
      dtot0<-symdynN(r = r, amu=amu, asd=asd, f=f, d=d,
                    d_sd=d_sd, d_cov=d_cov, N=max(Nlst),
                    sf=1, tmax=20, stochd = TRUE, stocht = TRUE, as.matrix = TRUE, amax = 0,
                    xstart = rnorm(max(Nlst), 0, sqrt(var_approx(r, f, d_sd))))
      x0<-unname(dtot0[nrow(dtot0),-c(1:2)])
      
      #simulation
      dtot<-symdynN(r = r, amu=amu, asd=asd, f=f, d=0,
                    d_sd=0, d_cov=0, N=max(Nlst),
                    sf=dtp, tmax=ntm, stochd = TRUE, stocht = TRUE, as.matrix = TRUE, amax = 0,
                    xstart = x0+rnorm(max(Nlst), 0, d_sd))
      m<-m+1
    }
    
    datout<-data.frame(dtot)
    rs<-abs(rowSums(datout[,-c(1:2)]))
    
    #raw r estimate
    rest<-(-log(rs[-1]/rs[1])/seq(dtp, ntm, by=dtp))
    rmat[ii,]<-rest
    
    if(ii/10 == floor(ii/10)) {
      cat(paste(round(ii/niter,2), ";"))
    }
  }
  write.csv(rmat, "datout/estmat_ec_r.csv", row.names = FALSE)
} else {
  rmat<-read.csv("datout/estmat_ec_r.csv")
}

padj<-c(-1, 0.02, 1.4)
qtlu<-0.975

pdf("figures/ecological_grain_time_effect.pdf", width=3, height=3, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(1,1), mar=c(3,4,1,1), oma=c(0.5,1.5,0,0))
  Asq<-seq(1, max(Nlst), by=0.1)
  Asq_small<-seq(1, max(Nlst))
  
  #r vs. sampling frequency
  rpuse<-seq(1, 100, by=3)
  tmlst<-c(seq(0.1, ntm, by=0.1)[rpuse])
  rmat2<-abs(rmat[,rpuse])
  tmp<-t(apply(rmat2, 2, function(x) quantile(x, pnorm(-1:1), na.rm=T)))
  pltqt(rep(tmlst, each=nrow(rmat2)), as.matrix(unlist(c(rmat2))), "", domod=FALSE, do_N = FALSE, plog = "", xlab = "", ylim=c(0,tmp[1,3]), jfac = 1, cluse = collst_attributes[2])
  lines(tmlst, tmp[,2], lwd=1.5, col=1, lty=2)
  abline(h=r, lty=3, lwd=1.5, col="black")
  mtext(expression(paste("resilience, ", italic(r))), 2, line=3.2)
  mtext(expression(paste("temporal scale")), 1, line=2.3)
  #title("a.", line=padj[1], adj=padj[2]+0.06, cex.main=padj[3])
dev.off()

pdf("figures/ecological_grain.pdf", width=3, height=6, colormodel = "cmyk", useDingbats = FALSE)
  par(mfcol=c(3,1), mar=c(2,4,1,1), oma=c(2,1.5,0,0))

  #r vs. number of species
  ps<-is.finite(estmat[,"r_naive"])
  #scaling reqlationship for variance
  var_est<-Asq*mean(estmat[which(estmat[,"N"]==max(Nlst)),"var_sp"],na.rm=T)+(Asq^2-Asq)*mean(estmat[which(estmat[,"N"]==max(Nlst)),"cov_sp"],na.rm=T)
  #scaling relationship for d_sd
  sd_N<-data.frame(Asq, dsd=sqrt((d_sd)^2*Asq*(1+(Asq-1)*d_cov/(d_sd)^2)))
  nps<-15 #which position to expand around
  #r at reference scale
  r_s1<-mean(estmat[,"r_naive"][estmat[,"N"]==nps])
  #variance at reference scale
  xmu_1<-sqrt(var_est[Asq==nps])
  
  #scaling relationship for r
  #mean species-level abundance and disturbance variance and covariance
  mu_varsp<-mean(estmat$var_sp[estmat$N==30],na.rm=T)
  mu_covsp<-mean(estmat$cov_sp[estmat$N==30],na.rm=T)
  mu_dsd<-mean(estmat$var_dist[estmat$N==30],na.rm=T)
  mu_covdist<-mean(estmat$cov_dist[estmat$N==30],na.rm=T) #basically zero
  
  #total summed covariance for abundance (V) and disturbance (D)
  Vsummed<-(mu_covsp*(Asq_small^2-Asq_small)+mu_varsp*Asq_small)
  Dsummed<-(mu_covdist*(Asq_small^2-Asq_small)+mu_dsd*Asq_small)

  #scaling relationship for r
  rest<-Dsummed/Vsummed/2
  ps<-which(is.finite(estmat[ps,"r_naive_median"]))
  
  pltqt(estmat[ps,"N"], estmat[ps,"r_naive_median"], "", domod=FALSE, do_N = FALSE, plog = "", xlab = "", ylim=c(0, 15), jfac = 1, cluse = collst_attributes[2])
  lines(Asq_small, rest, col=1, lwd=1.5, lty=2)
  
  mtext(expression(paste("resilience, ", italic(r))), 2, line=3.2)
  title("a.", line=padj[1], adj=padj[2], cex.main=padj[3])
  abline(h=r, lwd=1.5, lty=3, col="black")
  
  #d_sd
  #disturbance - follows same scaling relationship as for space
  pltqt(estmat[,"N"], estmat[,"d_sd_true"], "", sd_N, domod=FALSE, do_N = FALSE, plog = "", xlab = "", ylim=c(0, 1.8), jfac = 1, cluse = collst_attributes[1])
  lines(sd_N$Asq, sd_N$dsd, lwd=1.5, lty=3, col="black")
  
  title("b.", line=padj[1], adj=padj[2], cex.main=padj[3])
  mtext(expression(paste("resistance"^{-1}, ", ", italic(sigma))), 2, line=3.2)
  
  #var
  #scaling relationship for variance
  var_N<-data.frame(Asq, vr=var_approx(r,f,d_sd)*Asq*(1+(Asq-1)*d_cov/(d_sd)^2)) #null expectation (no interactions)
  pltqt(estmat[,"N"], estmat[,"var"], "", do_N = FALSE, domod=FALSE, plog = "y", mlog="", xlab = "", ylim=c(0.05,3.2), jfac = 1, cluse = collst_attributes[3])
  lines(var_N$Asq, var_N$vr, lwd=1.5, lty=3, col="black")
  title("c.", line=padj[1], adj=padj[2], cex.main=padj[3])
  mtext(expression(paste("variability, ", "var(", italic(x), ")")), 2, line=3.2)
  mtext(expression(paste("ecological scale")), 1, line = 2.5, outer=F, adj = 0.5)
  
  ps<-which(estmat[,"N"]==max(Nlst))
  lines(Asq, Asq*mean(estmat[ps,"var_sp"],na.rm=T)+(Asq^2-Asq)*mean(estmat[ps,"cov_sp"],na.rm=T), col=1, lwd=1.5, lty=2) #heuristic estimate
dev.off()
