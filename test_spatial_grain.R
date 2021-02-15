#plots Figure 4 and Figure S1 in the main text

#load functions and packages
source("functions.R")
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
d_cov<-d_cov0<-(d_sd)^2/2 #covariance in disturbance size among patches
rsd = 0.1 # standard deviation for r
Ksd = 0.1 # standard deviation for K

##################################
#simulate varying spatial grain, no dispersal
##################################
#set parameters
Alst<-1:30 #patch numbers to test
A<-c(2) #number of patches for dispersal test
Ilst<-seq(0, 2, by=0.1) #dispersal rates to test
niter<-1000 #number of iterations
Ifuse = 1 # dispersal rate

#only run if no saved simulation is available
if(sum(grep("estmat_sp_210209.csv", dir("datout/")))==0) {
  Alst<-1:30
  niter<-1000
  d_cov<-d_cov0
  estmat<-as.matrix(data.frame(iter=rep(niter, each=length(Alst)), tsmp=rep(Alst, niter),
                     N=NA, n=NA,
                     var=NA, f=NA, r_mean=NA, r_median=NA, d_sd_true = NA, meandist = NA))
  
  n<-1
  for(i in 1:niter) {
    #simulate all max(Alst) patches
    ruse = abs(rnorm(max(Alst), r, rsd))
    Kuse = abs(rnorm(max(Alst), K, Ksd))
    
    tmpout<-symdynN(r = ruse, amu=0, asd=0, f=f, d=d,
                  d_sd=d_sd, d_cov=d_cov, N=max(Alst),
                  sf=sf, tmax=tmax, stochd = TRUE, stocht = TRUE, as.matrix = TRUE, Ifrac=Ifuse, dffun = df_col, fullout = TRUE, Ksim = Kuse)
    dtot<-as.matrix(tmpout$datout)
    dtot<-dtot[dtot[,"time"]>20,] #throw out burn-in time
    
    dtot_mean = rowMeans(dtot[,-c(1, 2)])
    
    for(j in 1:length(Alst)) {
      #subsample a subset of patches
      smp = sample((1:max(Alst))+2, Alst[j])
      dtmp<-dtot[,c(1, 2, smp),drop=FALSE]
      estmat[n,"N"]<-Alst[j]
      estmat[n,"n"]<-nrow(dtmp)
      #observations are summed abundance across sampled patches
      xsum<-rowSums(dtmp[,-c(1:2),drop=FALSE])
      
      estmat[n,"meandist"]<-median(dtot_mean/(xsum/Alst[j]))
      
      #variance
      varest<-mean(xsum^2,na.rm=T)
      estmat[n,"var"]<-varest
      #mean disturbance frequencey
      estmat[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1])
  
      dt<-mean(diff(dtmp[,"time"]))
      
      # record true disturbance strength
      estmat[n,"d_sd_true"]<-sd(rowSums(tmpout$dquant[,smp-2,drop=FALSE])) #achieved disturbance standard deviations
      
      #estimate realized r
      disttm_forward<-which(c(dtmp[,"disturbed"][-1], 1)==0 & dtmp[,"disturbed"]==1)
      if(sum(disttm_forward)>0) {
        tmpr = suppressWarnings(log(xsum[disttm_forward+1]/xsum[disttm_forward])/dt)
        tmpr = tmpr[is.finite(tmpr)]
        estmat[n,"r_mean"]<-mean(tmpr)
        estmat[n,"r_median"]<-median(tmpr)
      }
      n<-n+1
    }
    if(i/10 == floor(i/10)) {
      cat(paste(round(i/niter,2), ";"))
    }
  }
  
  ##################################
  #simulate varying spatial grain, with varying dispersal
  ##################################
  set.seed(200209) #set random seed
  
  A<-c(2)
  Ilst<-seq(0, 2, by=0.2)
  d_cov<-0
  niter<-1000
  estmat_disp<-as.matrix(data.frame(iter=rep(niter, each=length(Ilst)), tsmp=rep(Ilst, niter),
                               n=NA, I=NA,
                               var=NA, f=NA, r_median=NA, d_sd_true=NA,
                               var_sp=NA, cor_sp=NA))
  n<-1
  for(i in 1:niter) {
    for(j in 1:length(Ilst)) {
      ruse = abs(rnorm(max(A), r, rsd))
      Kuse = K # abs(rnorm(max(A), K, Ksd))

      #simulate dynamics
      tmpout<-symdynN(r = ruse, amu=0, asd=0, f=f, d=d,
                    d_sd=d_sd, d_cov=d_cov, N=A,
                    sf=sf, tmax=tmax, stochd = TRUE, stocht = TRUE, as.matrix = TRUE, Ifrac=Ilst[j], dffun = df_col, fullout = TRUE, Ksim = Kuse)
      dtmp<-as.matrix(tmpout$datout)
      #remove first 20 steps for burn-in
      dtmp<-dtmp[dtmp[,"time"]>20,]
      
      estmat_disp[n,"I"]<-Ilst[j]
      estmat_disp[n,"n"]<-nrow(dtmp)
      
      #observation is summed abundance across the plots
      xsum<-rowSums(dtmp[,-c(1:2),drop=FALSE])
      
      dt<-mean(diff(dtmp[,"time"]))
      
      #variance
      varest<-mean(xsum^2,na.rm=T)
      estmat_disp[n,"var"]<-varest
      #average time between disturbances
      estmat_disp[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1])
      
      #patch-level variance and between-patch covariance
      estmat_disp[n,"var_sp"]<-mean(apply(dtmp[,-c(1:2),drop=FALSE], 2, function(x) mean(x^2)))
      cvtmp<-cov(dtmp[,-c(1:2),drop=FALSE])
      estmat_disp[n,"cor_sp"]<-mean(cvtmp[row(cvtmp)!=col(cvtmp)])
      
      # record true disturbance strength
      estmat_disp[n,"d_sd_true"]<-sd(rowSums(tmpout$dquant)) #achieved disturbance standard deviations
      
      #estimate realized r
      disttm_forward<-which(c(dtmp[,"disturbed"][-1], 1)==0 & dtmp[,"disturbed"]==1)
      if(sum(disttm_forward)>0) {
        tmpr = suppressWarnings(log(xsum[disttm_forward+1]/xsum[disttm_forward])/dt)
        tmpr = tmpr[is.finite(tmpr)]
        estmat_disp[n,"r_median"]<-median(tmpr)
      }
      n<-n+1
    }
    if(i/10 == floor(i/10)) {
      cat(paste(round(i/niter,2), ";"))
    }
  }
  
  write.csv(estmat, "datout/estmat_sp_210209.csv", row.names=F)
  write.csv(estmat_disp, "datout/estmat_sp_disp_210209.csv", row.names=F)
} else {
  estmat<-read.csv("datout/estmat_sp_210209.csv")
  estmat_disp<-read.csv("datout/estmat_sp_disp_210209.csv")
}

pdf("figures/spatial_grain_210209.pdf", width=3, height=6, colormodel = "cmyk", useDingbats = FALSE)
    par(mfcol=c(3,1), mar=c(2,4,1,1), oma=c(2,1.5,0,0))
    
    #r
    ps<-is.finite(estmat[,"r_median"])
    pltqt(estmat[ps,"N"], -estmat[ps,"r_median"], "", r, domod=FALSE,
          do_N = FALSE, plog = "", xlab = "spatial grain", ylim = c(0.89, 1.51),
          jfac=1, cluse = collst_attributes[2])
    title("a.", line=-1.2, adj=0.08, cex.main=1.4)
    mtext(expression(paste("resilience, ", italic(r))), 2, line=3.2)
    
    
    meandist_lst = tapply(estmat[,"meandist"], estmat[,"N"], function(x) median(x))
    meanr_lst = tapply(-estmat[,"r_median"], estmat[,"N"], function(x) mean(x))
    r_est = -meandist_lst*Ifuse+r+Ifuse
    
    lines(Alst, r_est, lty = 3)
    
    #d_sd
    Asq<-seq(1, max(Alst), length=1000)
    
    #scaling relationship for d_sd
    sd_N<-data.frame(A=Asq, d_sd=sqrt((d_sd)^2*Asq*(1+(Asq-1)*d_cov0/(d_sd)^2)))
    ps<-is.finite(estmat[,"d_sd_true"])
    #back-calculate from corrected r value
    pltqt(estmat[ps,"N"], estmat[ps,"d_sd_true"], "", sd_N, domod=FALSE, do_N = FALSE, plog = "", xlab = "spatial grain", ylim = c(0, 8), jfac=1, cluse = collst_attributes[1])
    title("b.", line=-1.2, adj=0.04, cex.main=1.4)
    mtext(expression(paste("resistance"^{-1}, ", ", italic(sigma))), 2, line=3.2)
    
    #variance
    #scaling relationship for variance
    var_N<-data.frame(Alst, var_approx(r_est,f,d_sd)*Alst*(1+(Alst-1)*d_cov0/(d_sd)^2))
    pltqt(estmat[,"N"], estmat[,"var"], "", truev = var_N, do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "spatial grain", ylim=c(0, 28), jfac=1, cluse = collst_attributes[3])
    title("c.", line=-1.2, adj=0.04, cex.main=1.4)
    mtext(expression(paste("invariance"^-1, ", ", "var(", italic(x), ")")), 2, line=3.2)
    
    mtext("spatial scale", 1, line = 2.5, outer=F, adj = 0.5)
dev.off()

pdf("figures/spatial_grain_dispersal_patch_level_210209.pdf", width=6, height=4, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,2), mar=c(2,4,1,1), oma=c(2,1,0,0))
  
  #r
  ps<-is.finite(estmat_disp[,"r_median"])
  pltqt(estmat_disp[ps,"I"], -estmat_disp[ps,"r_median"], "", r, domod=FALSE, do_N = FALSE,
        plog = "", xlab = "dispersal rate", ylim = c(0.9, 1.1), qtp = c(-1, 1),
        jfac = 1, cluse = collst_attributes[2])
  mtext(expression(paste("meta. resilience")), 2, line=3.1)
  title("a.", line=-1, adj=0.04, cex.main=1.2)

  
  #d_sd
  sd_N<-sqrt((d_sd)^2*A) #recall, d_cov=0
  ps<-is.finite(estmat_disp[,"d_sd_true"])
  pltqt(estmat_disp[ps,"I"], estmat_disp[ps,"d_sd_true"], "", sd_N, domod=FALSE, do_N = FALSE, plog = "", xlab = "dispersal rate", ylim=c(0.38, 0.52), jfac = 1, cluse = collst_attributes[1])
  #pltqt(estmat_disp[ps,"I"], sqrt(estmat_disp[ps,"r"]*(2*f*estmat_disp[ps,"var"])), "", sd_N, domod=FALSE, do_N = FALSE, plog = "", xlab = "dispersal rate", ylim=c(0.38, 0.52), jfac = 1, cluse = collst_attributes[1])
  mtext(expression(paste("meta. resistance"^{-1})), 2, line=3.1)
  title("b.", line=-1, adj=0.04, cex.main=1.2)
  
  #variance
  var_N<-var_approx(r,f,d_sd)*A #recall, d_cov=0
  pltqt(estmat_disp[,"I"], estmat_disp[,"var"], "", truev = var_N, do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "dispersal rate", ylim = c(0.07, 0.13), jfac = 1, cluse = collst_attributes[3])
  mtext(expression(paste("meta. invariance"^-1)), 2, line=3.1)
  title("c.", line=-1, adj=0.04, cex.main=1.2)
  
  
  pltqt(estmat_disp[,"I"], estmat_disp[,"var_sp"]*2, "", do_N = FALSE, domod=FALSE, plog = "", mlog="xy", xlab = "dispersal rate", modoffset = 1e-1, jfac = 1, ylim = c(-0.02, 0.14), cluse=collst_attributes[3])
  mtext(expression(paste("patch-level invariance"^-1)), 2, line=3.1)
  title("d.", line=-1, adj=0.02, cex.main=1.2)
  
  addqt(estmat_disp[,"I"], estmat_disp[,"cor_sp"]*(2^2-2), jfac = 1, cluse="darkblue", pltdens = NA, border=NA)
  abline(h=var_approx(r,f,d_sd), lwd=1, col="black", lty=3)
  
  text(1.5, 0.10, expression(paste(italic(A), paste("<var(", italic(x[i]), ")>"))), col=collst_attributes[3], cex=1.2)
  text(1.3, -0.03, expression(paste("(", italic(A)^2, "-", italic(A), ")", paste("<cov(", italic(x[i]), ",", italic(x[j]), ")>"))), col="darkblue", cex=1.2, pos = 3)
  
  mtext("dispersal rate", 1, line = 0.5, outer=TRUE)
dev.off()

